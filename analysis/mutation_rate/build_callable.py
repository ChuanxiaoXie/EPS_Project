#!/usr/bin/env python3
"""Build a merged callable BED from one Sentieon/GATK-style gVCF."""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import Counter
from pathlib import Path

try:
    from .common import load_config, thresholds
except ImportError:
    from common import load_config, thresholds


def integer_or_none(value: str) -> int | None:
    value = value.strip()
    return None if value in {"", "."} else int(float(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--gvcf", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.sample not in config["samples"]:
        raise ValueError(f"Sample is absent from configuration: {args.sample}")
    limits = thresholds(config)
    chromosomes = list(map(str, config["chromosomes"]))
    chromosome_set = set(chromosomes)
    out_bed = Path(args.out_bed)
    out_json = Path(args.out_json)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    query_format = r"%CHROM\t%POS\t%END[\t%GT\t%DP\t%GQ\t%MIN_DP]\n"
    command = [
        args.bcftools,
        "query",
        "-r",
        ",".join(chromosomes),
        "-s",
        args.sample,
        "-f",
        query_format,
        args.gvcf,
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True, encoding="utf-8")
    assert process.stdout is not None

    examined = accepted = raw_bp = merged_bp = 0
    per_chromosome_bp: Counter[str] = Counter()
    per_chromosome_records: Counter[str] = Counter()
    current_chromosome: str | None = None
    current_start = current_end = -1

    def flush(handle) -> None:
        nonlocal current_chromosome, current_start, current_end, merged_bp
        if current_chromosome is None:
            return
        handle.write(f"{current_chromosome}\t{current_start}\t{current_end}\n")
        length = current_end - current_start
        merged_bp += length
        per_chromosome_bp[current_chromosome] += length
        current_chromosome = None

    with out_bed.open("w", encoding="utf-8") as handle:
        for line in process.stdout:
            columns = [value.strip() for value in line.rstrip("\n").split("\t")]
            if len(columns) != 7:
                raise RuntimeError(f"Unexpected gVCF query row: {line[:200]}")
            chromosome, position_text, end_text, genotype, depth_text, gq_text, min_depth_text = columns
            if chromosome not in chromosome_set:
                continue
            examined += 1
            position = int(position_text)
            end = integer_or_none(end_text) or position
            depth = integer_or_none(depth_text)
            genotype_quality = integer_or_none(gq_text)
            minimum_depth = integer_or_none(min_depth_text)
            lower_depth = minimum_depth if minimum_depth is not None else depth
            if (
                genotype in {".", "./.", ".|."}
                or depth is None
                or lower_depth is None
                or genotype_quality is None
                or lower_depth < int(limits["callable_min_dp"])
                or depth > int(limits["callable_max_dp"])
                or genotype_quality < int(limits["callable_min_gq"])
            ):
                continue
            start = position - 1
            accepted += 1
            raw_bp += end - start
            per_chromosome_records[chromosome] += 1
            if current_chromosome == chromosome and start <= current_end:
                current_end = max(current_end, end)
            else:
                flush(handle)
                current_chromosome, current_start, current_end = chromosome, start, end
        flush(handle)

    return_code = process.wait()
    if return_code != 0:
        raise SystemExit(f"bcftools query failed with exit code {return_code}")
    payload = {
        "sample": args.sample,
        "input_gvcf": args.gvcf,
        "thresholds": limits,
        "chromosomes": chromosomes,
        "records_examined": examined,
        "records_accepted": accepted,
        "raw_callable_bp_before_merge": raw_bp,
        "merged_callable_bp": merged_bp,
        "per_chromosome_callable_bp": dict(per_chromosome_bp),
        "per_chromosome_accepted_records": dict(per_chromosome_records),
    }
    out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
