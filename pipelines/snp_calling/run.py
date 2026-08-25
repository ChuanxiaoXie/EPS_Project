#!/usr/bin/env python3
"""Render the complete per-sample and joint SNP-calling script set."""

from __future__ import annotations

import argparse
import csv
import shlex
import sys
from pathlib import Path, PurePosixPath


REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))

from eps_workflows.common import execute_script, load_json, require_keys, validate_identifier, write_script
from eps_workflows.sentieon_gvcf import render_sample
from eps_workflows.sentieon_joint_calling import render as render_joint


def generated_gvcf(output_root: str, sample_id: str) -> str:
    return str(PurePosixPath(output_root) / "gvcf" / sample_id / f"{sample_id}.g.vcf.gz")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gvcf-config", type=Path, required=True)
    parser.add_argument("--joint-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    gvcf_config = load_json(args.gvcf_config)
    joint_config = load_json(args.joint_config)
    require_keys(gvcf_config, ["output_root", "samples"], "gVCF configuration")
    samples = gvcf_config["samples"]
    if not isinstance(samples, list) or not samples:
        raise ValueError("samples must be a non-empty array")

    sample_dir = args.output_dir / "samples"
    sample_dir.mkdir(parents=True, exist_ok=True)
    scripts: list[Path] = []
    manifest_rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    for sample in samples:
        sample_id = validate_identifier(str(sample.get("sample_id", "")), "sample_id")
        if sample_id in seen:
            raise ValueError(f"Duplicate sample_id: {sample_id}")
        seen.add(sample_id)
        script = sample_dir / f"{sample_id}.sentieon_gvcf.sh"
        write_script(script, render_sample(gvcf_config, sample))
        scripts.append(script)
        manifest_rows.append((sample_id, generated_gvcf(str(gvcf_config["output_root"]), sample_id)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = args.output_dir / "generated.gvcfs.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "gvcf"])
        writer.writerows(manifest_rows)

    joint_script = args.output_dir / "sentieon_joint_calling.sh"
    write_script(joint_script, render_joint(joint_config, [row[1] for row in manifest_rows]))
    master = args.output_dir / "run_snp_calling_pipeline.sh"
    master_lines = [
        "#!/usr/bin/env bash",
        "set -Eeuo pipefail",
        "umask 002",
        'SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)',
        "",
    ]
    master_lines.extend(f'bash "$SCRIPT_DIR/samples/{script.name}"' for script in scripts)
    master_lines.extend(['bash "$SCRIPT_DIR/sentieon_joint_calling.sh"', ""])
    write_script(master, "\n".join(master_lines))
    print(master)
    if args.execute:
        execute_script(master)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
