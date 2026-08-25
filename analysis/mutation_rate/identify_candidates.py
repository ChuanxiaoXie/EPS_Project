#!/usr/bin/env python3
"""Identify first-appearance biallelic SNV candidates in a configured lineage."""

from __future__ import annotations

import argparse
import bisect
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    from .common import load_config, thresholds
except ImportError:
    from common import load_config, thresholds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--callable-dir", required=True)
    parser.add_argument("--mask-bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    return parser.parse_args()


def number_or_none(value: str) -> float | None:
    value = value.strip()
    return None if value in {"", "."} else float(value)


def parse_call(field: str) -> dict[str, Any]:
    parts = field.split(":")
    genotype = parts[0] if parts else "."
    allele_depth_text = parts[1] if len(parts) > 1 else "."
    depth = int(float(parts[2])) if len(parts) > 2 and parts[2] not in {"", "."} else None
    genotype_quality = int(float(parts[3])) if len(parts) > 3 and parts[3] not in {"", "."} else None
    if "." in genotype:
        code = None
    else:
        alleles = [int(value) for value in re.split(r"[/|]", genotype)]
        code = 0 if all(value == 0 for value in alleles) else (2 if all(value > 0 for value in alleles) else 1)
    allele_depths: list[int] = []
    if allele_depth_text not in {"", "."}:
        try:
            allele_depths = [int(value) for value in allele_depth_text.split(",")]
        except ValueError:
            allele_depths = []
    reference_depth = allele_depths[0] if allele_depths else None
    alternate_depth = allele_depths[1] if len(allele_depths) > 1 else None
    denominator = (reference_depth or 0) + (alternate_depth or 0)
    allele_balance = alternate_depth / denominator if alternate_depth is not None and denominator else None
    return {
        "gt": genotype,
        "code": code,
        "dp": depth,
        "gq": genotype_quality,
        "ref_ad": reference_depth,
        "alt_ad": alternate_depth,
        "ab": allele_balance,
    }


def read_bed(path: Path) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                chromosome, start, end, *_ = line.rstrip().split("\t")
                raw[chromosome].append((int(start), int(end)))
    return {
        chromosome: ([interval[0] for interval in sorted(intervals)], sorted(intervals))
        for chromosome, intervals in raw.items()
    }


def in_bed(index, chromosome: str, position: int) -> bool:
    if chromosome not in index:
        return False
    starts, intervals = index[chromosome]
    offset = bisect.bisect_right(starts, position - 1) - 1
    return offset >= 0 and intervals[offset][0] <= position - 1 < intervals[offset][1]


def keys_in_bed(path: Path, wanted: dict[str, list[int]]) -> set[tuple[str, int]]:
    """Stream a large sorted BED while retaining only queried positions."""
    found: set[tuple[str, int]] = set()
    offsets: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chromosome, start_text, end_text, *_ = line.rstrip().split("\t")
            positions = wanted.get(chromosome, [])
            if not positions:
                continue
            start, end = int(start_text), int(end_text)
            offset = offsets[chromosome]
            while offset < len(positions) and positions[offset] - 1 < start:
                offset += 1
            while offset < len(positions) and positions[offset] - 1 < end:
                found.add((chromosome, positions[offset]))
                offset += 1
            offsets[chromosome] = offset
    return found


def reference_high_confidence(call: dict[str, Any], limits: dict[str, float | int]) -> bool:
    return (
        call["code"] == 0
        and call["dp"] is not None
        and int(limits["callable_min_dp"]) <= int(call["dp"]) <= int(limits["callable_max_dp"])
        and call["gq"] is not None
        and int(call["gq"]) >= int(limits["callable_min_gq"])
    )


def target_strict(call: dict[str, Any], limits: dict[str, float | int]) -> bool:
    if (
        call["code"] not in {1, 2}
        or call["dp"] is None
        or not int(limits["callable_min_dp"]) <= int(call["dp"]) <= int(limits["callable_max_dp"])
        or call["gq"] is None
        or int(call["gq"]) < int(limits["callable_min_gq"])
        or call["alt_ad"] is None
        or int(call["alt_ad"]) < int(limits["target_min_alt_reads"])
        or call["ab"] is None
    ):
        return False
    balance = float(call["ab"])
    if call["code"] == 1:
        return float(limits["heterozygous_ab_min"]) <= balance <= float(limits["heterozygous_ab_max"])
    return balance >= float(limits["homozygous_alt_ab_min"])


def target_relaxed(call: dict[str, Any], limits: dict[str, float | int]) -> bool:
    return (
        call["code"] in {1, 2}
        and call["dp"] is not None
        and int(limits["callable_min_dp"]) <= int(call["dp"]) <= int(limits["callable_max_dp"])
        and call["alt_ad"] is not None
        and int(call["alt_ad"]) >= int(limits["target_min_alt_reads"])
    )


def transmission_consistent(
    calls: dict[str, dict[str, Any]], descendants: list[str], target: str, limits: dict[str, float | int], require_gq: bool
) -> bool:
    """Reject confident impossible regain or loss while ignoring uncertain descendants."""
    state = calls[target]["code"]
    for sample in descendants:
        call = calls[sample]
        high_confidence = (
            call["code"] is not None
            and call["dp"] is not None
            and int(limits["callable_min_dp"]) <= int(call["dp"]) <= int(limits["callable_max_dp"])
            and (
                not require_gq
                or (call["gq"] is not None and int(call["gq"]) >= int(limits["callable_min_gq"]))
            )
        )
        if not high_confidence:
            continue
        code = int(call["code"])
        if state == 0 and code != 0:
            return False
        if state == 2 and code != 2:
            return False
        if state == 1:
            state = code
    return True


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    limits = thresholds(config)
    samples = list(map(str, config["samples"]))
    chromosomes = list(map(str, config["chromosomes"]))
    accepted_filters = set(map(str, config.get("accepted_filters", ["PASS"])))
    allow_missing_mq = bool(config.get("allow_missing_site_mq", False))
    interval_by_target = {str(interval["target"]): interval for interval in config["intervals"]}
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mask = read_bed(Path(args.mask_bed))

    query_format = r"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum[\t%GT:%AD:%DP:%GQ]\n"
    command = [
        args.bcftools,
        "query",
        "-r",
        ",".join(chromosomes),
        "-s",
        ",".join(samples),
        "-f",
        query_format,
        args.vcf,
    ]
    process = subprocess.run(command, check=True, capture_output=True, text=True, encoding="utf-8")
    records: list[tuple[Any, ...]] = []
    positions: dict[str, list[int]] = defaultdict(list)
    for line in process.stdout.splitlines():
        columns = [value.strip() for value in line.split("\t")]
        if len(columns) != 11 + len(samples):
            raise RuntimeError(f"Unexpected VCF query column count: {len(columns)}")
        chromosome, position_text, reference, alternate, *rest = columns
        position = int(position_text)
        records.append((chromosome, position, reference, alternate, *rest))
        if len(reference) == 1 and len(alternate) == 1 and alternate not in {"*", "."}:
            positions[chromosome].append(position)
    for chromosome in positions:
        positions[chromosome].sort()

    callable_keys = {
        str(interval["name"]): keys_in_bed(
            Path(args.callable_dir) / f"{interval['name']}.callable.unmasked.bed", positions
        )
        for interval in config["intervals"]
    }
    rows: list[dict[str, Any]] = []
    for record in records:
        chromosome, position, reference, alternate = record[:4]
        quality, filter_value, qd_text, fs_text, mq_text, mq_rank_text, read_pos_text = record[4:11]
        sample_fields = record[11:]
        if len(reference) != 1 or len(alternate) != 1 or alternate in {"*", "."}:
            continue
        calls = {sample: parse_call(field) for sample, field in zip(samples, sample_fields)}
        matching_intervals: list[dict[str, Any]] = []
        for sample in samples:
            candidate_interval = interval_by_target.get(sample)
            if candidate_interval is None:
                continue
            ancestors = list(map(str, candidate_interval["ancestors"]))
            if calls[sample]["code"] in {1, 2} and all(calls[ancestor]["code"] == 0 for ancestor in ancestors):
                matching_intervals.append(candidate_interval)
        if not matching_intervals:
            continue
        interval = matching_intervals[0]
        parallel_shared = len({str(value.get("lineage", "")) for value in matching_intervals}) > 1
        interval_name = str(interval["name"])
        target = str(interval["target"])
        ancestors = list(map(str, interval["ancestors"]))
        descendants = list(map(str, interval.get("descendants", [])))
        if (chromosome, position) not in callable_keys[interval_name]:
            continue

        position_index = bisect.bisect_left(positions[chromosome], position)
        distances: list[int] = []
        if position_index > 0:
            distances.append(position - positions[chromosome][position_index - 1])
        if position_index + 1 < len(positions[chromosome]):
            distances.append(positions[chromosome][position_index + 1] - position)
        nearest_distance = min(distances) if distances else 10**12
        qd = number_or_none(qd_text)
        fs = number_or_none(fs_text)
        mq = number_or_none(mq_text)
        mq_rank = number_or_none(mq_rank_text)
        read_pos = number_or_none(read_pos_text)
        filter_pass = filter_value in accepted_filters
        site_strict = (
            filter_pass
            and qd is not None
            and qd >= float(limits["site_min_qd"])
            and fs is not None
            and fs <= float(limits["site_max_fs"])
            and (allow_missing_mq or (mq is not None and mq >= float(limits["site_min_mq"])))
            and (mq_rank is None or mq_rank >= float(limits["site_min_mq_rank_sum"]))
            and (read_pos is None or read_pos >= float(limits["site_min_read_pos_rank_sum"]))
            and nearest_distance >= int(limits["nearest_snv_distance"])
        )
        site_relaxed = (
            filter_pass
            and ((allow_missing_mq and mq is None) or (mq is not None and mq >= 40))
            and nearest_distance >= int(limits["nearest_snv_distance"])
        )
        ancestors_high = all(reference_high_confidence(calls[sample], limits) for sample in ancestors)
        strict_transmission = transmission_consistent(calls, descendants, target, limits, require_gq=True)
        relaxed_transmission = transmission_consistent(calls, descendants, target, limits, require_gq=False)
        strict_pre_bam = site_strict and ancestors_high and target_strict(calls[target], limits) and strict_transmission
        relaxed_pre_bam = site_relaxed and ancestors_high and target_relaxed(calls[target], limits) and relaxed_transmission
        if not (strict_pre_bam or relaxed_pre_bam):
            continue
        rows.append(
            {
                "candidate_id": f"{chromosome}:{position}:{reference}>{alternate}",
                "chrom": chromosome,
                "pos": position,
                "start0": position - 1,
                "end": position,
                "ref": reference,
                "alt": alternate,
                "interval": interval_name,
                "target_sample": target,
                "ancestor_samples": ",".join(ancestors),
                "parallel_branch_shared": "yes" if parallel_shared else "no",
                "qual": quality,
                "filter": filter_value,
                "QD": qd_text,
                "FS": fs_text,
                "MQ": mq_text,
                "MQRankSum": mq_rank_text,
                "ReadPosRankSum": read_pos_text,
                "nearest_biallelic_snv_distance": nearest_distance,
                "in_mask": "yes" if in_bed(mask, chromosome, position) else "no",
                "strict_pre_bam": "yes" if strict_pre_bam else "no",
                "relaxed_pre_bam": "yes" if relaxed_pre_bam else "no",
                "strict_transmission_consistent": "yes" if strict_transmission else "no",
                "relaxed_transmission_consistent": "yes" if relaxed_transmission else "no",
                "target_gt": calls[target]["gt"],
                "target_dp": calls[target]["dp"],
                "target_gq": calls[target]["gq"],
                "target_ref_ad": calls[target]["ref_ad"],
                "target_alt_ad": calls[target]["alt_ad"],
                "target_ab": calls[target]["ab"],
                "lineage_gt": ",".join(str(calls[sample]["gt"]) for sample in samples),
            }
        )

    fields = list(rows[0]) if rows else ["candidate_id"]
    with (outdir / "candidates.pre_bam.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (outdir / "candidates.pre_bam.bed").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['chrom']}\t{row['start0']}\t{row['end']}\t{row['candidate_id']}\n")
    print(f"pre_bam_candidates\t{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
