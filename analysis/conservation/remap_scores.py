#!/usr/bin/env python3
"""Remap position-based conservation scores through a coordinate map."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--position-map", required=True)
    parser.add_argument("--scores", required=True, help="Whitespace-delimited score file, optionally gzip-compressed")
    parser.add_argument("--output", required=True)
    parser.add_argument("--missing-value", default="NA")
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def load_scores(path: Path) -> dict[tuple[str, str], tuple[str, str]]:
    scores: dict[tuple[str, str], tuple[str, str]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Score row {line_number} has fewer than four columns")
            key = (fields[0], fields[1])
            if key in scores:
                raise ValueError(f"Duplicate score coordinate: {key[0]}:{key[1]}")
            scores[key] = (fields[2], fields[3])
    return scores


def parse_old_position(value: str) -> tuple[str, str]:
    fields = value.rsplit("_", 2)
    if len(fields) < 2:
        raise ValueError(f"Cannot parse mapped coordinate: {value}")
    return fields[0], fields[1]


def main() -> int:
    args = parse_args()
    scores = load_scores(Path(args.scores))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    found = missing = 0
    with Path(args.position_map).open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as target:
        target.write("chrom\tstart\tend\tscore1\tscore2\n")
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Map row {line_number} has fewer than four columns")
            key = parse_old_position(fields[3])
            values = scores.get(key)
            if values is None:
                values = (args.missing_value, args.missing_value)
                missing += 1
            else:
                found += 1
            target.write("\t".join([fields[0], fields[1], fields[2], values[0], values[1]]) + "\n")
    print(f"mapped\t{found}")
    print(f"missing\t{missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
