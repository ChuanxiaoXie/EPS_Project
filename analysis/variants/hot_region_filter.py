#!/usr/bin/env python3
"""Apply the manuscript hot-region SNP retention rule to a multi-sample VCF."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.variants.genotype_table import open_text
from analysis.variants.window_density import read_chromosome_sizes


Interval = tuple[int, int]
WindowRow = tuple[str, int, int, int, bool]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path, help="Input VCF or gzip-compressed VCF.")
    parser.add_argument("--chrom-sizes", required=True, type=Path, help="Two-column chromosome-size table.")
    parser.add_argument("--output-vcf", required=True, type=Path, help="Filtered VCF output; .gz uses gzip.")
    parser.add_argument("--hot-regions", required=True, type=Path, help="Merged hot-region BED output.")
    parser.add_argument("--window-table", required=True, type=Path, help="Auditable sliding-window count table.")
    parser.add_argument("--positions-output", required=True, type=Path, help="Retained CHROM/POS table.")
    parser.add_argument("--window-size", type=int, default=200_000)
    parser.add_argument("--step-size", type=int, default=20_000)
    parser.add_argument(
        "--hot-count-threshold",
        type=int,
        default=100,
        help="A window is hot when its SNP count is strictly greater than this value.",
    )
    parser.add_argument(
        "--accepted-hot-filter",
        action="append",
        help="FILTER value retained inside hot regions; default: PASS. Repeat for multiple values.",
    )
    parser.add_argument(
        "--legacy-pos-as-bed-start",
        action="store_true",
        help="Reproduce the legacy POS/POS point-BED convention instead of converting VCF POS to POS-1.",
    )
    return parser.parse_args()


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def is_snv(reference: str, alternate: str) -> bool:
    alleles = alternate.split(",")
    return len(reference) == 1 and bool(alleles) and all(len(value) == 1 and value not in {".", "*"} for value in alleles)


def read_snv_positions(
    path: Path, known_chromosomes: set[str], legacy_pos_as_bed_start: bool = False
) -> dict[str, list[int]]:
    positions: dict[str, list[int]] = defaultdict(list)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: malformed VCF row")
            chrom, raw_position, _identifier, reference, alternate = fields[:5]
            if chrom in known_chromosomes and is_snv(reference, alternate):
                position0 = int(raw_position) if legacy_pos_as_bed_start else int(raw_position) - 1
                if position0 < 0:
                    raise ValueError(f"{path}:{line_number}: VCF positions must be one-based")
                positions[chrom].append(position0)
    for values in positions.values():
        values.sort()
    return dict(positions)


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def select_hot_regions(
    positions: dict[str, list[int]],
    chromosome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    count_threshold: int,
) -> tuple[list[WindowRow], dict[str, list[Interval]]]:
    if window_size <= 0 or step_size <= 0:
        raise ValueError("Window and step sizes must be positive")
    if count_threshold < 0:
        raise ValueError("Hot-region count threshold cannot be negative")
    rows: list[WindowRow] = []
    hot_regions: dict[str, list[Interval]] = {}
    for chrom, chrom_length in chromosome_sizes.items():
        chrom_positions = positions.get(chrom, [])
        hot_windows: list[Interval] = []
        for start in range(0, chrom_length, step_size):
            end = min(start + window_size, chrom_length)
            count = bisect.bisect_left(chrom_positions, end) - bisect.bisect_left(chrom_positions, start)
            is_hot = count > count_threshold
            rows.append((chrom, start, end, count, is_hot))
            if is_hot:
                hot_windows.append((start, end))
            if end == chrom_length and start + step_size >= chrom_length:
                break
        hot_regions[chrom] = merge_intervals(hot_windows)
    return rows, hot_regions


def position_in_intervals(position0: int, intervals: list[Interval]) -> bool:
    starts = [start for start, _end in intervals]
    index = bisect.bisect_right(starts, position0) - 1
    return index >= 0 and position0 < intervals[index][1]


def write_windows(path: Path, rows: list[WindowRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "snp_count", "is_hot"])
        for chrom, start, end, count, is_hot in rows:
            writer.writerow([chrom, start, end, count, "yes" if is_hot else "no"])


def write_hot_regions(path: Path, regions: dict[str, list[Interval]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for chrom, intervals in regions.items():
            for start, end in intervals:
                handle.write(f"{chrom}\t{start}\t{end}\n")


def filter_vcf(
    input_path: Path,
    output_path: Path,
    positions_path: Path,
    hot_regions: dict[str, list[Interval]],
    accepted_hot_filters: set[str],
) -> dict[str, int]:
    counts = {"input_records": 0, "input_snvs": 0, "retained_snvs": 0, "hot_rejected_snvs": 0}
    positions_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(input_path) as source, open_output(output_path) as output, positions_path.open(
        "w", encoding="utf-8", newline=""
    ) as positions_handle:
        positions_writer = csv.writer(positions_handle, delimiter="\t", lineterminator="\n")
        positions_writer.writerow(["chrom", "pos"])
        for line_number, line in enumerate(source, start=1):
            if line.startswith("#"):
                output.write(line)
                continue
            if not line.strip():
                continue
            counts["input_records"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{input_path}:{line_number}: malformed VCF row")
            chrom, raw_position, _identifier, reference, alternate, _qual, filter_value = fields[:7]
            if not is_snv(reference, alternate):
                continue
            counts["input_snvs"] += 1
            position = int(raw_position)
            in_hot_region = position_in_intervals(position - 1, hot_regions.get(chrom, []))
            if in_hot_region and filter_value not in accepted_hot_filters:
                counts["hot_rejected_snvs"] += 1
                continue
            output.write(line if line.endswith("\n") else line + "\n")
            positions_writer.writerow([chrom, position])
            counts["retained_snvs"] += 1
    return counts


def main() -> int:
    args = parse_args()
    chromosome_sizes = read_chromosome_sizes(args.chrom_sizes)
    positions = read_snv_positions(args.vcf, set(chromosome_sizes), args.legacy_pos_as_bed_start)
    windows, hot_regions = select_hot_regions(
        positions,
        chromosome_sizes,
        args.window_size,
        args.step_size,
        args.hot_count_threshold,
    )
    write_windows(args.window_table, windows)
    write_hot_regions(args.hot_regions, hot_regions)
    counts = filter_vcf(
        args.vcf,
        args.output_vcf,
        args.positions_output,
        hot_regions,
        set(args.accepted_hot_filter or ["PASS"]),
    )
    summary = {
        **counts,
        "hot_region_bp": sum(end - start for values in hot_regions.values() for start, end in values),
        "hot_region_count": sum(len(values) for values in hot_regions.values()),
        "window_size": args.window_size,
        "step_size": args.step_size,
        "hot_count_rule": f"> {args.hot_count_threshold}",
        "density_coordinate_mode": "legacy_POS" if args.legacy_pos_as_bed_start else "VCF_POS_minus_1",
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
