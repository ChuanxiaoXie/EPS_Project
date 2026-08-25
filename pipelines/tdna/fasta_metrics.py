#!/usr/bin/env python3
"""Calculate basic contiguity metrics from a FASTA assembly."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_lengths(path: Path) -> list[int]:
    """Return sequence lengths while validating the minimal FASTA structure."""
    lengths: list[int] = []
    current_length: int | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                if current_length is not None:
                    lengths.append(current_length)
                current_length = 0
                continue
            if current_length is None:
                raise ValueError(
                    f"Sequence data precedes the first FASTA header at line {line_number}"
                )
            invalid = set(line.upper()).difference("ACGTURYSWKMBDHVN.-")
            if invalid:
                raise ValueError(
                    f"Unsupported FASTA symbols at line {line_number}: "
                    f"{''.join(sorted(invalid))}"
                )
            current_length += len(line)
    if current_length is not None:
        lengths.append(current_length)
    if not lengths:
        raise ValueError("The FASTA file contains no sequences")
    return lengths


def calculate_metrics(lengths: list[int]) -> dict[str, int]:
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    cumulative = 0
    n50 = 0
    for length in ordered:
        cumulative += length
        if cumulative * 2 >= total:
            n50 = length
            break
    return {
        "contigs": len(ordered),
        "total_bp": total,
        "N50": n50,
        "longest": ordered[0],
        "contigs_ge500": sum(length >= 500 for length in ordered),
        "contigs_ge1000": sum(length >= 1000 for length in ordered),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    metrics = calculate_metrics(read_lengths(args.fasta))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        for name, value in metrics.items():
            handle.write(f"{name}\t{value}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
