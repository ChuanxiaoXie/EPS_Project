#!/usr/bin/env python3
"""Render Sentieon joint genotyping and manuscript SNP hard filtering."""

from __future__ import annotations

import argparse
import csv
import shlex
from pathlib import Path
from typing import Any

from .common import execute_script, load_json, require_keys, write_script


def q(value: object) -> str:
    return shlex.quote(str(value))


def read_gvcf_manifest(path: Path) -> list[str]:
    """Read a two-column sample_id/gvcf manifest and reject ambiguous rows."""
    gvcfs: list[str] = []
    sample_ids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"sample_id", "gvcf"}.issubset(reader.fieldnames):
            raise ValueError("The gVCF manifest must contain sample_id and gvcf columns")
        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            gvcf = (row.get("gvcf") or "").strip()
            if not sample_id or not gvcf:
                raise ValueError(f"Incomplete gVCF manifest row at line {line_number}")
            if sample_id in sample_ids:
                raise ValueError(f"Duplicate sample_id in gVCF manifest: {sample_id}")
            sample_ids.add(sample_id)
            gvcfs.append(gvcf)
    if not gvcfs:
        raise ValueError("The gVCF manifest contains no samples")
    return gvcfs


def render(config: dict[str, Any], gvcfs: list[str]) -> str:
    """Return a non-destructive shell script for the confirmed publication flow."""
    require_keys(
        config,
        [
            "reference",
            "sentieon",
            "gatk",
            "license_server",
            "output_root",
            "joint_vcf_name",
        ],
        "configuration",
    )
    if not gvcfs:
        raise ValueError("At least one gVCF is required")

    threads = int(config.get("threads", 20))
    if threads < 1:
        raise ValueError("threads must be at least 1")
    output_root = str(config["output_root"]).rstrip("/")
    joint_vcf = f"{output_root}/{config['joint_vcf_name']}"
    snp_vcf = f"{output_root}/{config.get('snp_vcf_name', 'joint.snps.vcf.gz')}"
    filtered_vcf = f"{output_root}/{config.get('filtered_vcf_name', 'joint.snps.filtered.vcf.gz')}"
    variants = " \\\n+  ".join(f"-v {q(path)}" for path in gvcfs)

    # The thresholds reproduce the final manuscript source. They remain
    # configurable so a rerun can explicitly document any reviewed correction.
    thresholds = {
        "qd": float(config.get("min_qd", 2.0)),
        "fs": float(config.get("max_fs", 40.0)),
        "mq": float(config.get("min_mq", 40.0)),
        "mq_rank_sum": float(config.get("min_mq_rank_sum", -12.5)),
        "read_pos_rank_sum": float(config.get("min_read_pos_rank_sum", -4.0)),
        "haplotype_score": float(config.get("haplotype_score_threshold", 13.0)),
    }
    haplotype_operator = str(config.get("haplotype_score_operator", ">"))
    if haplotype_operator not in {"<", ">"}:
        raise ValueError("haplotype_score_operator must be '<' or '>'")
    expressions = {
        "qd": f"QD < {thresholds['qd']}",
        "fs": f"FS > {thresholds['fs']}",
        "mq": f"MQ < {thresholds['mq']}",
        "mq_rank_sum": f"MQRankSum < {thresholds['mq_rank_sum']}",
        "read_pos_rank_sum": f"ReadPosRankSum < {thresholds['read_pos_rank_sum']}",
        "haplotype_score": (
            f"haplotype_score {haplotype_operator} {thresholds['haplotype_score']}"
        ),
    }

    lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        "umask 002",
        "",
        "# Inputs and software locations are supplied by the public configuration.",
        f"export SENTIEON_LICENSE={q(config['license_server'])}",
        f"mkdir -p {q(output_root)}",
        "",
        "# Jointly genotype all per-sample gVCFs with the method used in the source analysis.",
        f"{q(config['sentieon'])} driver -r {q(config['reference'])} -t {threads} --algo GVCFtyper \\",
        f"  {variants} {q(joint_vcf)}",
        "",
        "# Retain SNPs before applying the manuscript hard-filter expressions.",
        f"{q(config['gatk'])} SelectVariants -R {q(config['reference'])} -V {q(joint_vcf)} \\",
        f"  --select-type-to-include SNP -O {q(snp_vcf)}",
        f"{q(config['gatk'])} VariantFiltration -V {q(snp_vcf)} -O {q(filtered_vcf)} \\",
        f"  --filter-name LowQD --filter-expression {q(expressions['qd'])} \\",
        f"  --filter-name HighFS --filter-expression {q(expressions['fs'])} \\",
        f"  --filter-name LowMQ --filter-expression {q(expressions['mq'])} \\",
        f"  --filter-name LowMQRankSum --filter-expression {q(expressions['mq_rank_sum'])} \\",
        f"  --filter-name LowReadPosRankSum --filter-expression {q(expressions['read_pos_rank_sum'])} \\",
        f"  --filter-name HaplotypeScoreFilter --filter-expression {q(expressions['haplotype_score'])}",
        "",
        f"printf '%s\\n' {q(filtered_vcf)}",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--gvcf-manifest", type=Path, required=True)
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    config = load_json(args.config)
    script_text = render(config, read_gvcf_manifest(args.gvcf_manifest))
    write_script(args.script, script_text)
    print(args.script)
    if args.execute:
        execute_script(args.script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
