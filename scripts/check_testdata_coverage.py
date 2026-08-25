#!/usr/bin/env python3
"""Verify that every public analysis or pipeline entry point has test fixtures."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
MATRIX = REPOSITORY / "testdata" / "test_matrix.tsv"
EXCLUDED_HELPERS = {
    "analysis/figures/export_utils.R",
    "analysis/mutation_rate/common.py",
}


def discover_entry_points() -> set[str]:
    entry_points: set[str] = set()
    for root in ("analysis", "pipelines", "workflows"):
        for path in (REPOSITORY / root).rglob("*"):
            if not path.is_file() or path.suffix not in {".py", ".R", ".sh"}:
                continue
            relative = path.relative_to(REPOSITORY).as_posix()
            if path.name == "__init__.py" or relative in EXCLUDED_HELPERS:
                continue
            entry_points.add(relative)
    return entry_points


def is_git_ignored(path: str) -> bool:
    if not (REPOSITORY / ".git").exists():
        return False
    result = subprocess.run(
        ["git", "check-ignore", "--quiet", "--", path],
        cwd=REPOSITORY,
        check=False,
    )
    return result.returncode == 0


def main() -> int:
    with MATRIX.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required_columns = {"entry_point", "mode", "fixtures", "note"}
    if not rows or not required_columns <= set(rows[0]):
        raise ValueError(f"{MATRIX} requires columns: {', '.join(sorted(required_columns))}")
    registered: dict[str, dict[str, str]] = {}
    errors: list[str] = []
    allowed_modes = {"executed", "rendered", "unit_tested", "optional_integration"}
    for row in rows:
        entry_point = row["entry_point"].strip()
        if entry_point in registered:
            errors.append(f"duplicate matrix entry: {entry_point}")
            continue
        registered[entry_point] = row
        if row["mode"] not in allowed_modes:
            errors.append(f"invalid mode for {entry_point}: {row['mode']}")
        for fixture in row["fixtures"].split(";"):
            fixture = fixture.strip()
            if fixture and fixture != "-" and not (REPOSITORY / fixture).is_file():
                errors.append(f"missing fixture for {entry_point}: {fixture}")
            elif fixture and fixture != "-" and is_git_ignored(fixture):
                errors.append(f"fixture is excluded from Git for {entry_point}: {fixture}")
    discovered = discover_entry_points()
    for entry_point in sorted(discovered - set(registered)):
        errors.append(f"unregistered public entry point: {entry_point}")
    for entry_point in sorted(set(registered) - discovered):
        errors.append(f"matrix references a missing entry point: {entry_point}")
    if errors:
        raise SystemExit("Test-data coverage check failed:\n- " + "\n- ".join(errors))
    counts = {mode: sum(row["mode"] == mode for row in rows) for mode in sorted(allowed_modes)}
    print(f"test-data coverage passed: {len(rows)} public entry points")
    print("; ".join(f"{mode}={count}" for mode, count in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
