#!/usr/bin/env python3
"""Render or execute the parameterized T-DNA SOAPdenovo workflow."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
from pathlib import Path


REQUIRED = ("project", "sample", "output_root", "bam_primary", "reference", "tdna_fasta")


def load_config(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    missing = [key for key in REQUIRED if not config.get(key)]
    if missing:
        raise ValueError(f"Missing required keys: {', '.join(missing)}")
    return config


def build_environment(config: dict[str, object]) -> dict[str, str]:
    mapping = {
        "PROJECT": config["project"],
        "SAMPLE": config["sample"],
        "SOAP_ROOT": config["output_root"],
        "BAM_PRIMARY": config["bam_primary"],
        "REF": config["reference"],
        "TDNA": config["tdna_fasta"],
        "SAMTOOLS": config.get("samtools", "samtools"),
        "CONDA_BIN": config.get("conda", "conda"),
        "FALLBACK_TOOL_BIN": config.get("fallback_tool_bin", ""),
        "TOOL_ENV": config.get("tool_env", f"{config['output_root']}/envs/assembly"),
        "NSLOTS": str(config.get("threads", 24)),
        "EXPECTED_TDNA_BP": config.get("expected_tdna_bp", ""),
        "TDNA_SCOPE_NOTE": config.get("tdna_scope_note", "Not specified"),
        "NORMALIZE_BAM_HEADER": "1" if config.get("normalize_bam_header", False) else "0",
    }
    if config.get("bam_secondary"):
        mapping["BAM_SECONDARY"] = config["bam_secondary"]
    if config.get("input_note"):
        mapping["INPUT_NOTE"] = config["input_note"]
    return {key: str(value) for key, value in mapping.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    environment = build_environment(config)
    pipeline = Path(__file__).with_name("pipeline.sh").resolve()
    rendered = " ".join(
        [*(f"{key}={shlex.quote(value)}" for key, value in environment.items()), "bash", shlex.quote(str(pipeline))]
    )
    print(rendered)
    if args.execute:
        process_environment = os.environ.copy()
        process_environment.update(environment)
        subprocess.run(["bash", str(pipeline)], env=process_environment, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
