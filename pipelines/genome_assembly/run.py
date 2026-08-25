#!/usr/bin/env python3
"""Render the recoverable genome scaffolding and assembly-assessment pipeline."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path, PurePosixPath
from typing import Any


REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))

from eps_workflows.common import execute_script, load_json, require_keys, write_script
from eps_workflows.hic_workflow import render as render_hic


def arguments(values: list[object]) -> str:
    return shlex.join(str(value) for value in values)


def render(config: dict[str, Any], hic_script: Path | None = None) -> str:
    require_keys(config, ["assessment_assembly", "output_root"], "configuration")
    output_root = str(config["output_root"])
    python = str(config.get("python", "python3"))
    minimum_length = int(config.get("minimum_sequence_length", 0))
    if minimum_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        "umask 002",
        'SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)',
        'REPOSITORY_ROOT=${REPOSITORY_ROOT:-.}',
        "",
        "# The supplied assessment assembly is an explicit pipeline input.",
        f"mkdir -p {shlex.quote(output_root)}",
    ]
    if hic_script is not None:
        lines.extend(
            [
                "",
                "# Optional Hi-C scaffolding stage.",
                f'bash "$SCRIPT_DIR/{hic_script.name}"',
            ]
        )

    lines.extend(
        [
            "",
            "# Contiguity and composition statistics.",
            arguments([python]) + ' "$REPOSITORY_ROOT/analysis/assembly/assembly_stats.py" ' + arguments(
                [
                    "--assembly",
                    config["assessment_assembly"],
                    "--output",
                    str(PurePosixPath(output_root) / "assembly_stats.tsv"),
                    "--minimum-length",
                    minimum_length,
                ]
            ),
        ]
    )

    merqury = config.get("merqury", {})
    if merqury.get("enabled", False):
        require_keys(merqury, ["reads", "merqury_root"], "merqury configuration")
        reads = merqury["reads"]
        if not isinstance(reads, list) or not reads:
            raise ValueError("merqury.reads must be a non-empty array")
        merqury_command: list[object] = []
        for read in reads:
            merqury_command.extend(["--read", read])
        merqury_command.extend(
            [
                "--assembly",
                config["assessment_assembly"],
                "--output",
                str(PurePosixPath(output_root) / "merqury"),
                "--merqury-root",
                merqury["merqury_root"],
                "--prefix",
                merqury.get("prefix", "assembly"),
            ]
        )
        if "kmer" in merqury:
            merqury_command.extend(["--kmer", int(merqury["kmer"])])
        lines.extend(
            [
                "",
                "# K-mer completeness and consensus quality.",
                'bash "$REPOSITORY_ROOT/workflows/assembly/run_merqury.sh" ' + arguments(merqury_command),
            ]
        )

    busco = config.get("busco", {})
    if busco.get("enabled", False):
        require_keys(busco, ["input", "lineage", "mode", "run_name"], "BUSCO configuration")
        busco_command: list[object] = [
            "--input",
            busco["input"],
            "--lineage",
            busco["lineage"],
            "--mode",
            busco["mode"],
            "--output",
            str(PurePosixPath(output_root) / "busco"),
            "--run-name",
            busco["run_name"],
            "--threads",
            int(busco.get("threads", 15)),
        ]
        if busco.get("offline", False):
            busco_command.append("--offline")
        lines.extend(
            [
                "",
                "# BUSCO completeness with explicit mode and lineage.",
                'bash "$REPOSITORY_ROOT/workflows/assembly/run_busco.sh" ' + arguments(busco_command),
            ]
        )

    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--script-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    config = load_json(args.config)
    args.script_dir.mkdir(parents=True, exist_ok=True)
    hic_script: Path | None = None
    if config.get("hic_config"):
        hic_config_path = Path(str(config["hic_config"]))
        if not hic_config_path.is_absolute():
            hic_config_path = args.config.parent / hic_config_path
        hic_config = load_json(hic_config_path)
        hic_script = args.script_dir / "hic_scaffolding.sh"
        write_script(hic_script, render_hic(hic_config))

    pipeline_script = args.script_dir / "run_genome_assembly_pipeline.sh"
    write_script(pipeline_script, render(config, hic_script))
    print(pipeline_script)
    if args.execute:
        execute_script(pipeline_script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
