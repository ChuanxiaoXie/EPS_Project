#!/usr/bin/env python3
"""Render a resumable, dependency-aware genome-annotation stage runner."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path
from typing import Any

from .common import execute_script, load_json, require_keys, validate_identifier, write_script


def topological_order(stages: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_name: dict[str, dict[str, Any]] = {}
    for stage in stages:
        require_keys(stage, ["name", "command"], "stage")
        name = validate_identifier(str(stage["name"]), "stage name")
        if name in by_name:
            raise ValueError(f"Duplicate stage name: {name}")
        by_name[name] = stage

    pending = set(by_name)
    completed: set[str] = set()
    ordered: list[dict[str, Any]] = []
    while pending:
        ready = sorted(
            name
            for name in pending
            if set(map(str, by_name[name].get("depends_on", []))) <= completed
        )
        if not ready:
            unresolved = {name: by_name[name].get("depends_on", []) for name in sorted(pending)}
            raise ValueError(f"Unknown or cyclic stage dependencies: {unresolved}")
        for name in ready:
            ordered.append(by_name[name])
            pending.remove(name)
            completed.add(name)
    return ordered


def render(config: dict[str, Any]) -> str:
    require_keys(config, ["work_dir", "stages"], "configuration")
    if not isinstance(config["stages"], list) or not config["stages"]:
        raise ValueError("stages must be a non-empty JSON array")
    ordered = topological_order(config["stages"])
    work_dir = str(config["work_dir"])
    state_dir = f"{work_dir}/.workflow_state"
    log_dir = f"{work_dir}/logs"

    lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        "umask 002",
        f"WORK_DIR={shlex.quote(work_dir)}",
        f"STATE_DIR={shlex.quote(state_dir)}",
        f"LOG_DIR={shlex.quote(log_dir)}",
        'mkdir -p "$WORK_DIR" "$STATE_DIR" "$LOG_DIR"',
        "",
        "run_stage() {",
        "  local name=$1 command=$2",
        '  if [[ -s "$STATE_DIR/$name.done" ]]; then echo "SKIP $name"; return; fi',
        '  printf "%s\\tSTART\\t%s\\n" "$(date --iso-8601=seconds)" "$name" >> "$STATE_DIR/progress.tsv"',
        '  bash -lc "$command" >"$LOG_DIR/$name.stdout.log" 2>"$LOG_DIR/$name.stderr.log"',
        '  printf "%s\\n" "$(date --iso-8601=seconds)" > "$STATE_DIR/$name.done"',
        '  printf "%s\\tDONE\\t%s\\n" "$(date --iso-8601=seconds)" "$name" >> "$STATE_DIR/progress.tsv"',
        "}",
        "",
    ]
    for stage in ordered:
        name = str(stage["name"])
        for input_path in stage.get("inputs", []):
            lines.append(f"test -s {shlex.quote(str(input_path))} || {{ echo 'Missing input for {name}: {input_path}' >&2; exit 2; }}")
        lines.append(f"run_stage {shlex.quote(name)} {shlex.quote(str(stage['command']))}")
        for output_path in stage.get("outputs", []):
            lines.append(f"test -s {shlex.quote(str(output_path))} || {{ echo 'Missing output from {name}: {output_path}' >&2; exit 3; }}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    write_script(args.script, render(load_json(args.config)))
    print(args.script)
    if args.execute:
        execute_script(args.script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
