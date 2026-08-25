#!/usr/bin/env python3
"""Summarize variant density in fixed, zero-filled genomic windows."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, type=Path, help="Tab-separated variant table with a header.")
    parser.add_argument("--chrom-sizes", required=True, type=Path, help="Two-column chromosome-size table.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path.")
    parser.add_argument("--chrom-column", default="chrom", help="Chromosome column in the variant table.")
    parser.add_argument("--position-column", default="pos", help="One-based position column in the variant table.")
    parser.add_argument("--group-column", help="Optional column used to calculate separate density tracks.")
    parser.add_argument("--window-size", type=int, default=100_000, help="Window size in base pairs.")
    return parser.parse_args()


def read_chromosome_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and length")
            chrom, raw_length = fields[:2]
            length = int(raw_length)
            if length <= 0:
                raise ValueError(f"{path}:{line_number}: chromosome length must be positive")
            if chrom in sizes:
                raise ValueError(f"{path}:{line_number}: duplicate chromosome {chrom!r}")
            sizes[chrom] = length
    if not sizes:
        raise ValueError(f"No chromosome sizes found in {path}")
    return sizes


def count_variants(
    path: Path,
    sizes: dict[str, int],
    chrom_column: str,
    position_column: str,
    group_column: str | None,
    window_size: int,
) -> tuple[Counter[tuple[str, int, str]], set[str]]:
    counts: Counter[tuple[str, int, str]] = Counter()
    groups: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {chrom_column, position_column}
        if group_column:
            required.add(group_column)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing variant-table columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            chrom = row[chrom_column]
            if chrom not in sizes:
                raise ValueError(f"{path}:{line_number}: chromosome {chrom!r} has no declared size")
            position = int(row[position_column])
            if position < 1 or position > sizes[chrom]:
                raise ValueError(f"{path}:{line_number}: position {position} is outside {chrom}")
            group = row[group_column] if group_column else "all"
            groups.add(group)
            window_index = (position - 1) // window_size
            counts[(chrom, window_index, group)] += 1
    if not groups:
        groups.add("all")
    return counts, groups


def write_density(
    path: Path,
    sizes: dict[str, int],
    counts: Counter[tuple[str, int, str]],
    groups: set[str],
    window_size: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "group", "variant_count", "variants_per_mb"])
        for chrom, chrom_length in sizes.items():
            for window_index in range(math.ceil(chrom_length / window_size)):
                start = window_index * window_size
                end = min(start + window_size, chrom_length)
                window_mb = (end - start) / 1_000_000
                for group in sorted(groups):
                    count = counts[(chrom, window_index, group)]
                    writer.writerow([chrom, start, end, group, count, f"{count / window_mb:.8g}"])


def main() -> None:
    args = parse_args()
    if args.window_size <= 0:
        raise ValueError("--window-size must be positive")
    sizes = read_chromosome_sizes(args.chrom_sizes)
    counts, groups = count_variants(
        args.variants,
        sizes,
        args.chrom_column,
        args.position_column,
        args.group_column,
        args.window_size,
    )
    write_density(args.output, sizes, counts, groups, args.window_size)


if __name__ == "__main__":
    main()
