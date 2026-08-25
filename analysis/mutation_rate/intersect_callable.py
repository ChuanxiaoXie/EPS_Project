#!/usr/bin/env python3
"""Intersect callable regions for configured pedigree intervals."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

try:
    from .common import bed_bp, load_config
except ImportError:
    from common import bed_bp, load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--callable-dir", required=True)
    parser.add_argument("--mask-bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bedtools", default="bedtools")
    return parser.parse_args()


def run_to_file(command: list[str], output: Path) -> None:
    with output.open("w", encoding="utf-8") as handle:
        subprocess.run(command, stdout=handle, check=True)


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    callable_dir = Path(args.callable_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for interval in config["intervals"]:
        name = str(interval["name"])
        required = list(map(str, interval["required_callable_samples"]))
        current = callable_dir / f"{required[0]}.callable.bed"
        if not current.is_file():
            raise FileNotFoundError(current)
        temporary: list[Path] = []
        final_unmasked = outdir / f"{name}.callable.unmasked.bed"
        if len(required) == 1:
            final_unmasked.write_bytes(current.read_bytes())
        for step, sample in enumerate(required[1:], start=1):
            sample_bed = callable_dir / f"{sample}.callable.bed"
            if not sample_bed.is_file():
                raise FileNotFoundError(sample_bed)
            is_last = step == len(required) - 1
            output = final_unmasked if is_last else outdir / f".{name}.step{step}.bed"
            run_to_file([args.bedtools, "intersect", "-a", str(current), "-b", str(sample_bed)], output)
            if current in temporary:
                current.unlink()
                temporary.remove(current)
            if not is_last:
                temporary.append(output)
            current = output
        masked = outdir / f"{name}.callable.masked.bed"
        run_to_file([args.bedtools, "subtract", "-a", str(final_unmasked), "-b", args.mask_bed], masked)
        unmasked_bp = bed_bp(final_unmasked)
        masked_bp = bed_bp(masked)
        generations = int(interval["generations"])
        row: dict[str, object] = {
                "interval": name,
                "target_sample": interval["target"],
                "generations": generations,
                "required_callable_samples": ",".join(required),
                "callable_bp_unmasked": unmasked_bp,
                "callable_bp_masked": masked_bp,
                "site_generations_unmasked": unmasked_bp * generations,
                "site_generations_masked": masked_bp * generations,
        }
        for optional_key in ("lineage", "parent_sample", "child_sample"):
            if optional_key in interval:
                row[optional_key] = interval[optional_key]
        rows.append(row)

    table = outdir / "callable_opportunities.tsv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "mask_bed": args.mask_bed,
        "intervals": rows,
        "total_site_generations_unmasked": sum(int(row["site_generations_unmasked"]) for row in rows),
        "total_site_generations_masked": sum(int(row["site_generations_masked"]) for row in rows),
    }
    (outdir / "callable_opportunities.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
