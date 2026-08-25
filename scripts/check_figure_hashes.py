#!/usr/bin/env python3
"""Verify that every curated binary figure matches the release manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


FIGURE_SUFFIXES = {".pdf", ".png"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("figures/SHA256SUMS"))
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = args.manifest if args.manifest.is_absolute() else root / args.manifest

    expected: dict[str, str] = {}
    for line_number, raw_line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise ValueError(f"Invalid checksum manifest row at line {line_number}")
        expected[fields[1].strip()] = fields[0].lower()

    actual_paths = {
        path.relative_to(root).as_posix()
        for path in (root / "figures").rglob("*")
        if path.is_file() and path.suffix.lower() in FIGURE_SUFFIXES
    }
    if actual_paths != set(expected):
        missing = sorted(set(expected) - actual_paths)
        unlisted = sorted(actual_paths - set(expected))
        raise ValueError(f"Figure manifest mismatch; missing={missing}, unlisted={unlisted}")

    failures = []
    for relative_path, expected_hash in sorted(expected.items()):
        observed = sha256(root / relative_path)
        if observed != expected_hash:
            failures.append(f"{relative_path}: expected {expected_hash}, observed {observed}")
    if failures:
        print("Figure checksum verification failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print(f"Verified {len(expected)} curated figure files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
