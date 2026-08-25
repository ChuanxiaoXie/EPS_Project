"""Shared configuration and interval helpers for mutation-rate analysis."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
DEFAULT_THRESHOLDS: dict[str, float | int] = {
    "callable_min_dp": 5,
    "callable_max_dp": 40,
    "callable_min_gq": 30,
    "site_min_qd": 2.0,
    "site_max_fs": 40.0,
    "site_min_mq": 50.0,
    "site_min_mq_rank_sum": -12.5,
    "site_min_read_pos_rank_sum": -4.0,
    "nearest_snv_distance": 5,
    "target_min_alt_reads": 5,
    "heterozygous_ab_min": 0.25,
    "heterozygous_ab_max": 0.75,
    "homozygous_alt_ab_min": 0.90,
    "ancestor_min_ref_reads": 5,
    "ancestor_max_alt_reads": 1,
    "ancestor_max_ab": 0.05,
    "max_other_base_fraction": 0.05,
    "bam_min_mq": 50,
    "bam_min_bq": 30,
}


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        config = json.load(handle)
    samples = config.get("samples")
    intervals = config.get("intervals")
    chromosomes = config.get("chromosomes")
    if not isinstance(samples, list) or not samples:
        raise ValueError("Configuration requires a non-empty samples list")
    if len(samples) != len(set(samples)):
        raise ValueError("Sample identifiers must be unique")
    if any(not isinstance(sample, str) or not IDENTIFIER.fullmatch(sample) for sample in samples):
        raise ValueError("Sample identifiers must be filesystem-safe")
    if not isinstance(chromosomes, list) or not chromosomes:
        raise ValueError("Configuration requires a non-empty chromosomes list")
    if not isinstance(intervals, list) or not intervals:
        raise ValueError("Configuration requires a non-empty intervals list")
    names: set[str] = set()
    for interval in intervals:
        if not isinstance(interval, dict):
            raise ValueError("Each interval must be an object")
        required = {"name", "target", "generations", "ancestors", "required_callable_samples"}
        missing = sorted(required - set(interval))
        if missing:
            raise ValueError(f"Interval is missing keys: {', '.join(missing)}")
        name = str(interval["name"])
        if not IDENTIFIER.fullmatch(name) or name in names:
            raise ValueError(f"Invalid or duplicate interval name: {name}")
        names.add(name)
        referenced = [
            interval["target"],
            *interval["ancestors"],
            *interval.get("descendants", []),
            *interval["required_callable_samples"],
        ]
        unknown = sorted(set(map(str, referenced)) - set(samples))
        if unknown:
            raise ValueError(f"Interval {name} references unknown samples: {', '.join(unknown)}")
        if int(interval["generations"]) < 1:
            raise ValueError(f"Interval {name} must span at least one generation")
    return config


def thresholds(config: dict[str, Any]) -> dict[str, float | int]:
    values = dict(DEFAULT_THRESHOLDS)
    supplied = config.get("thresholds", {})
    if not isinstance(supplied, dict):
        raise ValueError("thresholds must be an object")
    unknown = sorted(set(supplied) - set(values))
    if unknown:
        raise ValueError(f"Unknown thresholds: {', '.join(unknown)}")
    values.update(supplied)
    return values


def bed_bp(path: str | Path) -> int:
    total = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3:
                raise ValueError(f"Invalid BED row in {path}: {line.rstrip()}")
            total += int(columns[2]) - int(columns[1])
    return total
