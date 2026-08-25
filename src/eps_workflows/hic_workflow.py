#!/usr/bin/env python3
"""Render a non-destructive Juicer and 3D-DNA workflow."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path, PurePosixPath
from typing import Any

from .common import execute_script, load_json, require_keys, validate_identifier, write_script


def q(value: object) -> str:
    return shlex.quote(str(value))


def join(root: str, *parts: str) -> str:
    return str(PurePosixPath(root).joinpath(*parts))


def render(config: dict[str, Any]) -> str:
    require_keys(
        config,
        [
            "sample_id",
            "reference",
            "read1",
            "read2",
            "output_root",
            "bwa",
            "python",
            "juicer_root",
            "three_d_dna_root",
        ],
        "configuration",
    )
    sample = validate_identifier(str(config["sample_id"]), "sample_id")
    output = str(config["output_root"])
    reference_link = join(output, "reference", "genome.fa")
    read1_link = join(output, "fastq", f"{sample}_R1.fastq.gz")
    read2_link = join(output, "fastq", f"{sample}_R2.fastq.gz")
    restriction = str(config.get("restriction_site", "DpnII"))
    site_file = join(output, "restriction", f"{sample}_{restriction}.txt")
    sizes_file = join(output, "restriction", f"{sample}.chrom.{restriction}.sizes")
    threads = int(config.get("threads", 60))
    juicer_script = join(str(config["juicer_root"]), "scripts", "juicer_v1.sh")
    site_generator = join(str(config["juicer_root"]), "misc", "generate_site_positions.py")
    three_d = join(str(config["three_d_dna_root"]), "run-asm-pipeline.sh")
    merged = join(output, "aligned", "merged_nodups.txt")
    scaffold_dir = join(output, "3d-dna")

    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -Eeuo pipefail",
            "umask 002",
            "",
            "# Existing regular files are never overwritten by link creation.",
            "safe_link() {",
            "  local source=$1 destination=$2",
            "  if [[ -L \"$destination\" ]]; then ln -sfn -- \"$source\" \"$destination\";",
            "  elif [[ -e \"$destination\" ]]; then echo \"Refusing to replace $destination\" >&2; return 2;",
            "  else ln -s -- \"$source\" \"$destination\"; fi",
            "}",
            "",
            f"mkdir -p {q(join(output, 'fastq'))} {q(join(output, 'restriction'))} {q(join(output, 'reference'))} {q(scaffold_dir)}",
            f"safe_link {q(config['reference'])} {q(reference_link)}",
            f"safe_link {q(config['read1'])} {q(read1_link)}",
            f"safe_link {q(config['read2'])} {q(read2_link)}",
            f"{q(config['bwa'])} index {q(reference_link)}",
            f"cd {q(output)}",
            f"{q(config['python'])} {q(site_generator)} {q(restriction)} {q(sample)} {q(reference_link)}",
            f"mv -- {q(sample + '_' + restriction + '.txt')} {q(site_file)}",
            f"awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$NF}}' {q(site_file)} > {q(sizes_file)}",
            "",
            "# The workflow does not delete aligned/; reruns must use a new output root or an explicit backup.",
            f"{q(juicer_script)} -z {q(reference_link)} -p {q(sizes_file)} -y {q(site_file)} -d {q(output)} -D {q(config['juicer_root'])} -s {q(restriction)} -g {q(sample)} -t {threads} > {q(join(output, 'juicer.log'))} 2>&1",
            f"test -s {q(merged)}",
            f"cd {q(scaffold_dir)}",
            f"{q(three_d)} -r {int(config.get('review_rounds', 2))} {q(reference_link)} {q(merged)} > {q(join(scaffold_dir, '3d-dna.log'))} 2>&1",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    script = args.script
    write_script(script, render(load_json(args.config)))
    print(script)
    if args.execute:
        execute_script(script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
