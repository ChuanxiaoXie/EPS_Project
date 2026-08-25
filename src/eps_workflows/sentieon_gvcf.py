#!/usr/bin/env python3
"""Render one reproducible Sentieon alignment and gVCF script per sample."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path, PurePosixPath
from typing import Any

from .common import execute_script, load_json, require_keys, validate_identifier, write_script


def q(value: object) -> str:
    return shlex.quote(str(value))


def posix_join(root: str, *parts: str) -> str:
    return str(PurePosixPath(root).joinpath(*parts))


def render_sample(config: dict[str, Any], sample: dict[str, Any]) -> str:
    require_keys(sample, ["sample_id", "read1", "read2"], "sample")
    sample_id = validate_identifier(str(sample["sample_id"]), "sample_id")
    sentieon = str(config["sentieon"])
    reference = str(config["reference"])
    output_root = str(config["output_root"])
    threads = int(config.get("threads", 30))
    ploidy = int(sample.get("ploidy", config.get("ploidy", 2)))
    cleanup = bool(config.get("cleanup_sorted_bam", False))
    filter_bam = bool(config.get("filter_bam", False))
    minimum_mapping_quality = int(config.get("minimum_mapping_quality", 20))
    maximum_edit_distance = int(config.get("maximum_edit_distance", 1))
    maximum_cigar_operations = int(config.get("maximum_cigar_operations", 2))
    if min(threads, ploidy, minimum_mapping_quality, maximum_edit_distance, maximum_cigar_operations) < 0:
        raise ValueError("Thread, ploidy and BAM-filter settings must be non-negative")
    if threads < 1 or ploidy < 1 or maximum_cigar_operations < 1:
        raise ValueError("threads, ploidy and maximum_cigar_operations must be at least 1")
    samtools = str(config.get("samtools", "samtools"))

    alignment_dir = posix_join(output_root, "alignment", sample_id)
    gvcf_dir = posix_join(output_root, "gvcf", sample_id)
    sorted_bam = posix_join(alignment_dir, f"{sample_id}.sort.bam")
    dedup_bam = posix_join(alignment_dir, f"{sample_id}.rmdup.bam")
    filtered_bam = posix_join(alignment_dir, f"{sample_id}.rmdup.filtered.bam")
    score = posix_join(alignment_dir, f"{sample_id}.SCORE.gz")
    metrics = posix_join(alignment_dir, sample_id)
    gvcf = posix_join(gvcf_dir, f"{sample_id}.g.vcf.gz")
    rg = f"@RG\\tID:{sample_id}\\tLB:{sample.get('library', sample_id)}\\tSM:{sample_id}"

    lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        "umask 002",
        "",
        "# Every path and sample attribute comes from the JSON configuration.",
        f"export SENTIEON_LICENSE={q(config['license_server'])}",
        f"mkdir -p {q(alignment_dir)} {q(gvcf_dir)}",
        "",
        f"{q(sentieon)} bwa mem -k 32 -M -R {q(rg)} {q(reference)} {q(sample['read1'])} {q(sample['read2'])} -t {threads} | \\",
        f"  {q(sentieon)} util sort -r {q(reference)} -o {q(sorted_bam)} -t {threads} --sam2bam -i -",
        "",
        f"{q(sentieon)} driver -t {threads} -r {q(reference)} -i {q(sorted_bam)} \\",
        f"  --algo GCBias --summary {q(metrics + '.GC_SUMMARY.txt')} {q(metrics + '.GC_METRIC.txt')} \\",
        f"  --algo MeanQualityByCycle {q(metrics + '.MQ_METRIC.txt')} \\",
        f"  --algo QualDistribution {q(metrics + '.QD_METRIC.txt')} \\",
        f"  --algo InsertSizeMetricAlgo {q(metrics + '.IS_METRIC.txt')} \\",
        f"  --algo AlignmentStat {q(metrics + '.ALN_METRIC.txt')}",
        f"{q(sentieon)} plot QualDistribution -o {q(metrics + '.QD_METRIC.pdf')} {q(metrics + '.QD_METRIC.txt')}",
        f"{q(sentieon)} plot InsertSizeMetricAlgo -o {q(metrics + '.IS_METRIC.pdf')} {q(metrics + '.IS_METRIC.txt')}",
        "",
        f"{q(sentieon)} driver -t {threads} -r {q(reference)} -i {q(sorted_bam)} --algo LocusCollector --fun score_info {q(score)}",
        f"{q(sentieon)} driver -t {threads} -i {q(sorted_bam)} --algo Dedup --rmdup --score_info {q(score)} \\",
        f"  --metrics {q(metrics + '.DEDUP_METRIC.txt')} {q(dedup_bam)}",
    ]
    haplotyper_bam = dedup_bam
    if filter_bam:
        awk_program = (
            "BEGIN { OFS=\"\\t\" } "
            "/^@/ { print; next } "
            "{ nm=-1; for (i=12; i<=NF; i++) { "
            "if ($i ~ /^NM:i:/) { split($i, tag, \":\"); nm=tag[3]; break } } "
            "cigar=$6; operation_count=gsub(/[MIDNSHP=X]/, \"\", cigar); "
            "if (nm >= 0 && nm <= max_nm && operation_count <= max_ops) print }"
        )
        lines.extend(
            [
                "",
                "# Retain properly paired primary alignments with the manuscript BAM quality rules.",
                f"{q(samtools)} view -@ {threads} -h -F 4 -F 256 -q {minimum_mapping_quality} -f 2 -F 2048 {q(dedup_bam)} | \\",
                f"  awk -v max_nm={maximum_edit_distance} -v max_ops={maximum_cigar_operations} {q(awk_program)} | \\",
                f"  {q(samtools)} view -@ {threads} -b -o {q(filtered_bam)} -",
                f"{q(samtools)} index -@ {threads} {q(filtered_bam)}",
            ]
        )
        haplotyper_bam = filtered_bam
    lines.extend(
        [
            "",
            f"{q(sentieon)} driver -r {q(reference)} -t {threads} -i {q(haplotyper_bam)} \\",
            f"  --algo Haplotyper --ploidy {ploidy} --emit_conf=30 --call_conf=30 --emit_mode gvcf {q(gvcf)}",
        ]
    )
    if cleanup:
        lines.extend(
            [
                "",
                "# Cleanup is opt-in and only targets the intermediate created above.",
                f"rm -f -- {q(sorted_bam)} {q(sorted_bam + '.bai')}",
            ]
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--script-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    config = load_json(args.config)
    require_keys(config, ["reference", "sentieon", "license_server", "output_root", "samples"], "configuration")
    if not isinstance(config["samples"], list) or not config["samples"]:
        raise ValueError("samples must be a non-empty JSON array")

    seen: set[str] = set()
    scripts: list[Path] = []
    for sample in config["samples"]:
        if not isinstance(sample, dict):
            raise ValueError("Each sample entry must be a JSON object")
        sample_id = validate_identifier(str(sample.get("sample_id", "")), "sample_id")
        if sample_id in seen:
            raise ValueError(f"Duplicate sample_id: {sample_id}")
        seen.add(sample_id)
        script = args.script_dir / f"{sample_id}.sentieon_gvcf.sh"
        write_script(script, render_sample(config, sample))
        scripts.append(script)

    for script in scripts:
        print(script)
        if args.execute:
            execute_script(script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
