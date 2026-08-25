#!/usr/bin/env python3
"""Compare a normalized query VCF with a normalized simulation truth VCF."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import TextIO


Variant = tuple[str, int, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--query-vcf", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--false-positive-vcf-keys")
    parser.add_argument("--false-negative-vcf-keys")
    parser.add_argument("--accepted-filter", action="append")
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_variants(path: Path, accepted_filters: set[str]) -> set[Variant]:
    variants: set[Variant] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                raise ValueError(f"VCF row {line_number} has fewer than seven columns")
            chromosome, position, _identifier, reference, alternates, _quality, filter_value = fields[:7]
            if filter_value not in accepted_filters:
                continue
            for alternate in alternates.split(","):
                if alternate not in {".", "*"}:
                    variants.add((chromosome, int(position), reference.upper(), alternate.upper()))
    return variants


def variant_type(variant: Variant) -> str:
    return "SNV" if len(variant[2]) == len(variant[3]) == 1 else "INDEL_OR_COMPLEX"


def write_keys(path: str | None, variants: set[Variant]) -> None:
    if not path:
        return
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("chrom\tpos\tref\talt\n")
        for variant in sorted(variants, key=lambda row: (row[0], row[1], row[2], row[3])):
            handle.write("\t".join(map(str, variant)) + "\n")


def main() -> int:
    args = parse_args()
    accepted = set(args.accepted_filter or ["PASS", "."])
    truth = read_variants(Path(args.truth_vcf), accepted)
    query = read_variants(Path(args.query_vcf), accepted)
    true_positive = truth & query
    false_positive = query - truth
    false_negative = truth - query
    precision = len(true_positive) / len(query) if query else 0.0
    recall = len(true_positive) / len(truth) if truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    payload = {
        "normalization_requirement": "Both inputs must be left-normalized and split into biallelic records against the same reference.",
        "accepted_filters": sorted(accepted),
        "truth_variants": len(truth),
        "query_variants": len(query),
        "true_positive": len(true_positive),
        "false_positive": len(false_positive),
        "false_negative": len(false_negative),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "true_positive_by_type": dict(Counter(map(variant_type, true_positive))),
        "false_positive_by_type": dict(Counter(map(variant_type, false_positive))),
        "false_negative_by_type": dict(Counter(map(variant_type, false_negative))),
    }
    output = Path(args.output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_keys(args.false_positive_vcf_keys, false_positive)
    write_keys(args.false_negative_vcf_keys, false_negative)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
