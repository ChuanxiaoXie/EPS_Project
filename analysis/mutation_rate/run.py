#!/usr/bin/env python3
"""Render or execute the complete callable-genome mutation-rate workflow."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
from pathlib import Path

try:
    from .common import load_config, thresholds
except ImportError:
    from common import load_config, thresholds


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--gvcf-manifest", required=True, help="TSV columns: sample, gvcf")
    parser.add_argument("--bam-manifest", required=True, help="TSV columns: sample, bam")
    parser.add_argument("--joint-vcf", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--mask-bed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def read_manifest(path: str | Path, value_column: str) -> dict[str, str]:
    with Path(path).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", value_column}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"Manifest requires columns: sample, {value_column}")
    result: dict[str, str] = {}
    for row in rows:
        sample = row["sample"].strip()
        value = row[value_column].strip()
        if not sample or not value or sample in result:
            raise ValueError(f"Invalid or duplicate manifest sample: {sample}")
        result[sample] = value
    return result


def command_text(command: list[str]) -> str:
    return " ".join(shlex.quote(value) for value in command)


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    limits = thresholds(config)
    samples = list(map(str, config["samples"]))
    gvcfs = read_manifest(args.gvcf_manifest, "gvcf")
    bams = read_manifest(args.bam_manifest, "bam")
    for label, manifest in (("gVCF", gvcfs), ("BAM", bams)):
        missing = sorted(set(samples) - set(manifest))
        extra = sorted(set(manifest) - set(samples))
        if missing or extra:
            raise ValueError(f"{label} manifest mismatch; missing={missing}, extra={extra}")

    output = Path(args.output).resolve()
    callable_dir = output / "callable" / "samples"
    interval_dir = output / "callable" / "intervals"
    candidate_dir = output / "candidates"
    rate_dir = output / "mutation_rates"
    commands: list[tuple[list[str], Path | None]] = []
    for sample in samples:
        commands.append(
            (
                [
                    sys.executable,
                    str(SCRIPT_DIR / "build_callable.py"),
                    "--config",
                    args.config,
                    "--gvcf",
                    gvcfs[sample],
                    "--sample",
                    sample,
                    "--out-bed",
                    str(callable_dir / f"{sample}.callable.bed"),
                    "--out-json",
                    str(callable_dir / f"{sample}.callable.json"),
                    "--bcftools",
                    args.bcftools,
                ],
                None,
            )
        )
    commands.append(
        (
            [
                sys.executable,
                str(SCRIPT_DIR / "intersect_callable.py"),
                "--config",
                args.config,
                "--callable-dir",
                str(callable_dir),
                "--mask-bed",
                args.mask_bed,
                "--outdir",
                str(interval_dir),
                "--bedtools",
                args.bedtools,
            ],
            None,
        )
    )
    commands.append(
        (
            [
                sys.executable,
                str(SCRIPT_DIR / "identify_candidates.py"),
                "--config",
                args.config,
                "--vcf",
                args.joint_vcf,
                "--callable-dir",
                str(interval_dir),
                "--mask-bed",
                args.mask_bed,
                "--outdir",
                str(candidate_dir),
                "--bcftools",
                args.bcftools,
            ],
            None,
        )
    )
    mpileup = candidate_dir / "candidates.mpileup"
    commands.append(
        (
            [
                args.samtools,
                "mpileup",
                "-q",
                str(limits["bam_min_mq"]),
                "-Q",
                str(limits["bam_min_bq"]),
                "--ff",
                "3844",
                "--rf",
                "2",
                "-f",
                args.reference,
                "-l",
                str(candidate_dir / "candidates.pre_bam.bed"),
                *[bams[sample] for sample in samples],
            ],
            mpileup,
        )
    )
    validated = candidate_dir / "candidates.bam_validated.tsv"
    commands.append(
        (
            [
                sys.executable,
                str(SCRIPT_DIR / "parse_mpileup.py"),
                "--config",
                args.config,
                "--candidates",
                str(candidate_dir / "candidates.pre_bam.tsv"),
                "--mpileup",
                str(mpileup),
                "--output",
                str(validated),
            ],
            None,
        )
    )
    commands.append(
        (
            [
                sys.executable,
                str(SCRIPT_DIR / "calculate_rate.py"),
                "--opportunities",
                str(interval_dir / "callable_opportunities.tsv"),
                "--candidates",
                str(validated),
                "--output-tsv",
                str(rate_dir / "mutation_rates.tsv"),
                "--output-json",
                str(rate_dir / "mutation_rates.json"),
                "--output-lineage-tsv",
                str(rate_dir / "lineage_rates.tsv"),
            ],
            None,
        )
    )

    plan = [
        {"command": command_text(command), "stdout": str(stdout_path) if stdout_path else None}
        for command, stdout_path in commands
    ]
    print(json.dumps(plan, indent=2))
    if not args.execute:
        return 0
    for directory in (callable_dir, interval_dir, candidate_dir, rate_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for command, stdout_path in commands:
        if stdout_path is None:
            subprocess.run(command, check=True)
        else:
            with stdout_path.open("w", encoding="utf-8") as handle:
                subprocess.run(command, stdout=handle, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
