#!/usr/bin/env python3
"""Filter BLAST tabular alignments by identity and query coverage."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths", required=True, help="TSV columns: query_id, query_length")
    parser.add_argument("--blast", required=True, help="BLAST outfmt 6 table")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-query-coverage", type=float, default=80.0)
    return parser.parse_args()


def load_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Length row {line_number} has fewer than two columns")
            lengths[fields[0]] = int(fields[1])
    return lengths


def main() -> int:
    args = parse_args()
    lengths = load_lengths(Path(args.lengths))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    selected_queries: set[str] = set()
    selected_rows = 0
    with Path(args.blast).open(encoding="utf-8") as source, output.open("w", encoding="utf-8") as target:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 12:
                raise ValueError(f"BLAST row {line_number} has fewer than twelve columns")
            query = fields[0]
            query_length = lengths.get(query)
            if not query_length:
                continue
            identity = float(fields[2])
            aligned_length = int(fields[3])
            coverage = aligned_length / query_length * 100
            if identity >= args.min_identity and coverage >= args.min_query_coverage:
                target.write(line if line.endswith("\n") else line + "\n")
                selected_queries.add(query)
                selected_rows += 1
    print(f"selected_queries\t{len(selected_queries)}")
    print(f"selected_alignments\t{selected_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
