#!/usr/bin/env python3
"""Parse multi-sample samtools mpileup and validate pedigree candidates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

try:
    from .common import load_config, thresholds
except ImportError:
    from common import load_config, thresholds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--mpileup", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def count_bases(bases: str, alternate: str) -> dict[str, int]:
    counts = {"ref_fwd": 0, "ref_rev": 0, "alt_fwd": 0, "alt_rev": 0, "other": 0}
    offset = 0
    while offset < len(bases):
        symbol = bases[offset]
        if symbol == "^":
            offset += 2
            continue
        if symbol == "$":
            offset += 1
            continue
        if symbol in "+-":
            offset += 1
            length_start = offset
            while offset < len(bases) and bases[offset].isdigit():
                offset += 1
            insertion_length = int(bases[length_start:offset]) if offset > length_start else 0
            offset += insertion_length
            continue
        if symbol == ".":
            counts["ref_fwd"] += 1
        elif symbol == ",":
            counts["ref_rev"] += 1
        elif symbol == alternate.upper():
            counts["alt_fwd"] += 1
        elif symbol == alternate.lower():
            counts["alt_rev"] += 1
        elif symbol in "ACGTNacgtn*<>":
            counts["other"] += 1
        offset += 1
    return counts


def empty_counts() -> dict[str, int | float]:
    return {
        "reported_depth": 0,
        "ref": 0,
        "alt": 0,
        "ref_fwd": 0,
        "ref_rev": 0,
        "alt_fwd": 0,
        "alt_rev": 0,
        "other": 0,
        "ab": 0.0,
    }


def pileup_state(counts: dict[str, Any], limits: dict[str, float | int]) -> int | None:
    if (
        int(counts["ref"]) >= int(limits["ancestor_min_ref_reads"])
        and int(counts["alt"]) <= int(limits["ancestor_max_alt_reads"])
        and float(counts["ab"]) <= float(limits["ancestor_max_ab"])
    ):
        return 0
    if (
        int(counts["ref"]) >= int(limits["target_min_alt_reads"])
        and int(counts["alt"]) >= int(limits["target_min_alt_reads"])
        and float(limits["heterozygous_ab_min"]) <= float(counts["ab"]) <= float(limits["heterozygous_ab_max"])
    ):
        return 1
    if (
        int(counts["alt"]) >= int(limits["target_min_alt_reads"])
        and float(counts["ab"]) >= float(limits["homozygous_alt_ab_min"])
    ):
        return 2
    return None


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    limits = thresholds(config)
    samples = list(map(str, config["samples"]))
    intervals = {str(interval["name"]): interval for interval in config["intervals"]}
    with Path(args.candidates).open(encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    by_key = {(row["chrom"], int(row["pos"])): row for row in candidates}
    pileups: dict[tuple[str, int], dict[str, dict[str, Any]]] = {}
    with Path(args.mpileup).open(encoding="utf-8") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3 + 3 * len(samples):
                continue
            key = (columns[0], int(columns[1]))
            if key not in by_key:
                continue
            alternate = by_key[key]["alt"]
            sample_counts: dict[str, dict[str, Any]] = {}
            for sample_index, sample in enumerate(samples):
                depth = int(columns[3 + sample_index * 3])
                bases = columns[4 + sample_index * 3]
                counts: dict[str, Any] = count_bases(bases, alternate)
                counts["reported_depth"] = depth
                counts["ref"] = counts["ref_fwd"] + counts["ref_rev"]
                counts["alt"] = counts["alt_fwd"] + counts["alt_rev"]
                denominator = counts["ref"] + counts["alt"]
                counts["ab"] = counts["alt"] / denominator if denominator else 0.0
                sample_counts[sample] = counts
            pileups[key] = sample_counts

    output_rows: list[dict[str, Any]] = []
    for row in candidates:
        key = (row["chrom"], int(row["pos"]))
        counts = pileups.get(key, {sample: empty_counts() for sample in samples})
        target = row["target_sample"]
        ancestors = row["ancestor_samples"].split(",") if row.get("ancestor_samples") else []
        descendants = list(map(str, intervals[row["interval"]].get("descendants", [])))
        target_counts = counts[target]
        target_genotype = row["target_gt"].replace("|", "/")
        if target_genotype in {"0/1", "1/0"}:
            target_pass = (
                int(target_counts["alt"]) >= int(limits["target_min_alt_reads"])
                and int(target_counts["ref"]) >= int(limits["target_min_alt_reads"])
                and float(limits["heterozygous_ab_min"])
                <= float(target_counts["ab"])
                <= float(limits["heterozygous_ab_max"])
            )
        else:
            target_pass = (
                int(target_counts["alt"]) >= int(limits["target_min_alt_reads"])
                and float(target_counts["ab"]) >= float(limits["homozygous_alt_ab_min"])
            )
        maximum_other = max(
            1,
            int(float(limits["max_other_base_fraction"]) * max(1, int(target_counts["reported_depth"]))),
        )
        other_bases_pass = int(target_counts["other"]) <= maximum_other
        ancestors_pass = all(pileup_state(counts[sample], limits) == 0 for sample in ancestors)
        lineage_pass = True
        state = pileup_state(target_counts, limits)
        for sample in descendants:
            code = pileup_state(counts[sample], limits)
            if code is None:
                continue
            if state == 0 and code != 0:
                lineage_pass = False
            elif state == 2 and code != 2:
                lineage_pass = False
            elif state == 1:
                state = code
        bam_strict = target_pass and other_bases_pass and ancestors_pass and lineage_pass
        strand_pass = int(target_counts["alt_fwd"]) >= 1 and int(target_counts["alt_rev"]) >= 1
        bam_relaxed = (
            int(target_counts["alt"]) >= int(limits["target_min_alt_reads"])
            and all(int(counts[sample]["ref"]) >= int(limits["ancestor_min_ref_reads"]) for sample in ancestors)
        )
        output = dict(row)
        output.update(
            {
                "bam_strict_pass": "yes" if bam_strict else "no",
                "bam_strict_strand_pass": "yes" if bam_strict and strand_pass else "no",
                "bam_relaxed_pass": "yes" if bam_relaxed else "no",
                "bam_lineage_consistent": "yes" if lineage_pass else "no",
            }
        )
        for sample in samples:
            sample_counts = counts[sample]
            output[f"{sample}_pileup_dp"] = sample_counts["reported_depth"]
            output[f"{sample}_ref"] = sample_counts["ref"]
            output[f"{sample}_alt"] = sample_counts["alt"]
            output[f"{sample}_ab"] = f"{float(sample_counts['ab']):.6f}"
            output[f"{sample}_alt_fwd"] = sample_counts["alt_fwd"]
            output[f"{sample}_alt_rev"] = sample_counts["alt_rev"]
            output[f"{sample}_other"] = sample_counts["other"]
        output_rows.append(output)

    fields = list(output_rows[0]) if output_rows else ["candidate_id"]
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"candidates\t{len(output_rows)}")
    print(f"bam_strict_pass\t{sum(row['bam_strict_pass'] == 'yes' for row in output_rows)}")
    print(f"bam_strict_strand_pass\t{sum(row['bam_strict_strand_pass'] == 'yes' for row in output_rows)}")
    print(f"bam_relaxed_pass\t{sum(row['bam_relaxed_pass'] == 'yes' for row in output_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
