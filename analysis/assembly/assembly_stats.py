#!/usr/bin/env python3
"""Calculate reproducible contiguity and composition metrics from an assembly FASTA."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


VALID_BASES = frozenset("ACGTURYSWKMBDHVN.-")


@dataclass(frozen=True)
class SequenceSummary:
    length: int
    gc_bases: int
    canonical_bases: int
    n_bases: int


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def read_fasta(path: Path, minimum_length: int = 0) -> list[SequenceSummary]:
    """Stream a FASTA file and retain per-sequence summaries only."""
    if minimum_length < 0:
        raise ValueError("minimum_length cannot be negative")
    summaries: list[SequenceSummary] = []
    length = gc_bases = canonical_bases = n_bases = 0
    in_record = False

    def finish_record() -> None:
        if in_record and length >= minimum_length:
            summaries.append(SequenceSummary(length, gc_bases, canonical_bases, n_bases))

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                finish_record()
                length = gc_bases = canonical_bases = n_bases = 0
                in_record = True
                continue
            if not in_record:
                raise ValueError(f"Sequence data precedes the first FASTA header at line {line_number}")
            sequence = line.upper()
            invalid = set(sequence).difference(VALID_BASES)
            if invalid:
                raise ValueError(
                    f"Unsupported FASTA symbols at line {line_number}: {''.join(sorted(invalid))}"
                )
            length += len(sequence)
            gc_bases += sequence.count("G") + sequence.count("C")
            canonical_bases += sum(sequence.count(base) for base in "ACGT")
            n_bases += sequence.count("N")
    finish_record()
    if not summaries:
        raise ValueError("The FASTA file contains no sequences meeting the minimum length")
    return summaries


def nx_lx(lengths: list[int], fraction: float) -> tuple[int, int]:
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be within (0, 1]")
    ordered = sorted(lengths, reverse=True)
    target = sum(ordered) * fraction
    cumulative = 0
    for index, length in enumerate(ordered, start=1):
        cumulative += length
        if cumulative >= target:
            return length, index
    raise AssertionError("Nx calculation did not reach its target")


def calculate_metrics(summaries: list[SequenceSummary]) -> dict[str, int | float]:
    lengths = [summary.length for summary in summaries]
    n50, l50 = nx_lx(lengths, 0.5)
    n90, l90 = nx_lx(lengths, 0.9)
    canonical = sum(summary.canonical_bases for summary in summaries)
    gc_bases = sum(summary.gc_bases for summary in summaries)
    return {
        "sequence_count": len(lengths),
        "total_bp": sum(lengths),
        "longest_bp": max(lengths),
        "shortest_bp": min(lengths),
        "N50_bp": n50,
        "L50": l50,
        "N90_bp": n90,
        "L90": l90,
        "N_bases": sum(summary.n_bases for summary in summaries),
        "GC_percent_of_ACGT": 100.0 * gc_bases / canonical if canonical else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-length", type=int, default=0)
    args = parser.parse_args()

    metrics = calculate_metrics(read_fasta(args.assembly, args.minimum_length))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        for name, value in metrics.items():
            rendered = f"{value:.6f}" if isinstance(value, float) else str(value)
            handle.write(f"{name}\t{rendered}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
