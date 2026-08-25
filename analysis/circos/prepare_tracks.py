#!/usr/bin/env python3
"""Prepare the five 200-kb tracks used by the manuscript chromosome Circos plot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.variants.genotype_table import open_text
from analysis.variants.hot_region_filter import is_snv, merge_intervals, position_in_intervals
from analysis.variants.window_density import read_chromosome_sizes


Interval = tuple[int, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path, help="Selected SNP VCF used for the Circos tracks.")
    parser.add_argument("--genes-bed", required=True, type=Path, help="BED-like gene or transcript intervals.")
    parser.add_argument("--reference", required=True, type=Path, help="Reference FASTA used for G+C counts.")
    parser.add_argument("--chrom-sizes", required=True, type=Path, help="Two-column chromosome-size table.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--window-size", type=int, default=200_000)
    parser.add_argument("--karyotype-color", default="chr1", help="Circos color assigned to every chromosome.")
    parser.add_argument(
        "--legacy-pos-as-bed-start",
        action="store_true",
        help="Reproduce the legacy POS/POS point-BED convention instead of converting VCF POS to POS-1.",
    )
    return parser.parse_args()


def read_variant_positions(
    path: Path, chromosomes: set[str], legacy_pos_as_bed_start: bool = False
) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    all_snps: dict[str, list[int]] = defaultdict(list)
    gc_to_at: dict[str, list[int]] = defaultdict(list)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: malformed VCF row")
            chrom, raw_position, _identifier, reference, alternate = fields[:5]
            if chrom not in chromosomes or not is_snv(reference, alternate):
                continue
            position0 = int(raw_position) if legacy_pos_as_bed_start else int(raw_position) - 1
            if position0 < 0:
                raise ValueError(f"{path}:{line_number}: VCF positions must be one-based")
            all_snps[chrom].append(position0)
            if reference.upper() in {"G", "C"} and any(value.upper() in {"A", "T"} for value in alternate.split(",")):
                gc_to_at[chrom].append(position0)
    for collection in (all_snps, gc_to_at):
        for positions in collection.values():
            positions.sort()
    return dict(all_snps), dict(gc_to_at)


def read_gene_intervals(path: Path, chromosome_sizes: dict[str, int]) -> dict[str, list[Interval]]:
    genes: dict[str, list[Interval]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected chrom, start and end")
            chrom = fields[0]
            if chrom not in chromosome_sizes:
                continue
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start or end > chromosome_sizes[chrom]:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            genes[chrom].append((start, end))
    for intervals in genes.values():
        intervals.sort()
    return dict(genes)


def empty_track(chromosome_sizes: dict[str, int], window_size: int) -> dict[str, list[int]]:
    return {
        chrom: [0] * math.ceil(chrom_length / window_size)
        for chrom, chrom_length in chromosome_sizes.items()
    }


def count_points(
    chromosome_sizes: dict[str, int], positions: dict[str, list[int]], window_size: int
) -> dict[str, list[int]]:
    counts = empty_track(chromosome_sizes, window_size)
    for chrom, values in positions.items():
        for position0 in values:
            if position0 >= chromosome_sizes[chrom]:
                raise ValueError(f"Variant position {position0 + 1} exceeds declared length for {chrom}")
            counts[chrom][position0 // window_size] += 1
    return counts


def count_overlapping_features(
    chromosome_sizes: dict[str, int], genes: dict[str, list[Interval]], window_size: int
) -> dict[str, list[int]]:
    counts = empty_track(chromosome_sizes, window_size)
    for chrom, intervals in genes.items():
        for start, end in intervals:
            first_window = start // window_size
            last_window = (end - 1) // window_size
            for window_index in range(first_window, last_window + 1):
                counts[chrom][window_index] += 1
    return counts


def select_gene_snps(
    positions: dict[str, list[int]], genes: dict[str, list[Interval]]
) -> dict[str, list[int]]:
    selected: dict[str, list[int]] = {}
    for chrom, values in positions.items():
        merged_genes = merge_intervals(genes.get(chrom, []))
        selected[chrom] = [position for position in values if position_in_intervals(position, merged_genes)]
    return selected


def count_reference_gc(
    fasta_path: Path, chromosome_sizes: dict[str, int], window_size: int
) -> dict[str, list[int]]:
    counts = empty_track(chromosome_sizes, window_size)
    observed_lengths: dict[str, int] = {}
    current_chrom: str | None = None
    current_position = 0
    with open_text(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_chrom in chromosome_sizes:
                    observed_lengths[current_chrom] = current_position
                current_chrom = line[1:].strip().split()[0]
                current_position = 0
                continue
            sequence = line.strip().upper()
            if not sequence or current_chrom not in chromosome_sizes:
                continue
            offset = 0
            while offset < len(sequence):
                window_index = current_position // window_size
                if window_index >= len(counts[current_chrom]):
                    raise ValueError(f"Reference sequence {current_chrom} exceeds its declared length")
                take = min(len(sequence) - offset, (window_index + 1) * window_size - current_position)
                segment = sequence[offset : offset + take]
                counts[current_chrom][window_index] += segment.count("G") + segment.count("C")
                current_position += take
                offset += take
    if current_chrom in chromosome_sizes:
        observed_lengths[current_chrom] = current_position
    missing = set(chromosome_sizes).difference(observed_lengths)
    if missing:
        raise ValueError(f"Reference FASTA is missing chromosomes: {', '.join(sorted(missing))}")
    mismatched = {
        chrom: (observed_lengths[chrom], length)
        for chrom, length in chromosome_sizes.items()
        if observed_lengths[chrom] != length
    }
    if mismatched:
        details = ", ".join(f"{chrom}={observed}/{declared}" for chrom, (observed, declared) in mismatched.items())
        raise ValueError(f"Reference and chromosome-size lengths differ: {details}")
    return counts


def write_track(
    path: Path, chromosome_sizes: dict[str, int], values: dict[str, list[int]], window_size: int
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, chrom_length in chromosome_sizes.items():
            for window_index, value in enumerate(values[chrom]):
                start = window_index * window_size
                writer.writerow([chrom, start, min(start + window_size, chrom_length), value])


def write_karyotype(path: Path, chromosome_sizes: dict[str, int], color: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for chrom, chrom_length in chromosome_sizes.items():
            handle.write(f"chr\t-\t{chrom}\t{chrom}\t0\t{chrom_length}\t{color}\n")


def prepare_tracks(
    vcf_path: Path,
    genes_path: Path,
    fasta_path: Path,
    sizes_path: Path,
    output_dir: Path,
    window_size: int,
    karyotype_color: str = "chr1",
    legacy_pos_as_bed_start: bool = False,
) -> dict[str, int]:
    if window_size <= 0:
        raise ValueError("Window size must be positive")
    chromosome_sizes = read_chromosome_sizes(sizes_path)
    all_snps, gc_to_at = read_variant_positions(
        vcf_path, set(chromosome_sizes), legacy_pos_as_bed_start
    )
    genes = read_gene_intervals(genes_path, chromosome_sizes)
    tracks = {
        "1_wai.txt": count_points(chromosome_sizes, all_snps, window_size),
        "1_nei.txt": count_points(chromosome_sizes, select_gene_snps(all_snps, genes), window_size),
        "2_genedensity.txt": count_overlapping_features(chromosome_sizes, genes, window_size),
        "3_GC_wai.txt": count_reference_gc(fasta_path, chromosome_sizes, window_size),
        "3_GC2AT_nei.txt": count_points(chromosome_sizes, gc_to_at, window_size),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, values in tracks.items():
        write_track(output_dir / name, chromosome_sizes, values, window_size)
    write_karyotype(output_dir / "karyotype.txt", chromosome_sizes, karyotype_color)
    summary = {
        "chromosome_count": len(chromosome_sizes),
        "selected_snp_count": sum(len(values) for values in all_snps.values()),
        "gene_snp_count": sum(len(values) for values in select_gene_snps(all_snps, genes).values()),
        "gene_feature_count": sum(len(values) for values in genes.values()),
        "gc_to_at_snp_count": sum(len(values) for values in gc_to_at.values()),
        "window_size": window_size,
        "variant_coordinate_mode": "legacy_POS" if legacy_pos_as_bed_start else "VCF_POS_minus_1",
    }
    (output_dir / "track_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    args = parse_args()
    summary = prepare_tracks(
        args.vcf,
        args.genes_bed,
        args.reference,
        args.chrom_sizes,
        args.output_dir,
        args.window_size,
        args.karyotype_color,
        args.legacy_pos_as_bed_start,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
