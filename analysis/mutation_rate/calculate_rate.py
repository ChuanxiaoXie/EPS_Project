#!/usr/bin/env python3
"""Calculate mutation-rate scenarios and exact Poisson confidence intervals."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opportunities", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-lineage-tsv")
    return parser.parse_args()


def poisson_interval(count: int, alpha: float = 0.05) -> tuple[float, float]:
    try:
        from scipy.stats import chi2
    except ImportError as error:
        raise RuntimeError("scipy is required for exact Poisson confidence intervals") from error
    lower = 0.0 if count == 0 else 0.5 * float(chi2.ppf(alpha / 2, 2 * count))
    upper = 0.5 * float(chi2.ppf(1 - alpha / 2, 2 * (count + 1)))
    return lower, upper


def exact_rate_ratio(count_a: int, exposure_a: int, count_b: int, exposure_b: int) -> dict[str, object]:
    try:
        from scipy.stats import beta, binomtest
    except ImportError as error:
        raise RuntimeError("scipy is required for exact lineage comparisons") from error
    if exposure_a <= 0 or exposure_b <= 0:
        return {"rate_ratio": None, "exact95_low": None, "exact95_high": None, "p_value": None}
    total_count = count_a + count_b
    if total_count == 0:
        return {"rate_ratio": None, "exact95_low": 0.0, "exact95_high": None, "p_value": 1.0}
    null_probability = exposure_a / (exposure_a + exposure_b)
    p_value = float(binomtest(count_a, total_count, null_probability).pvalue)
    lower_probability = 0.0 if count_a == 0 else float(beta.ppf(0.025, count_a, total_count - count_a + 1))
    upper_probability = 1.0 if count_a == total_count else float(beta.ppf(0.975, count_a + 1, total_count - count_a))

    def convert(probability: float) -> float | str:
        return "Infinity" if probability >= 1 else probability / (1 - probability) * exposure_b / exposure_a

    if count_b == 0:
        rate_ratio: float | str | None = "Infinity" if count_a else None
    else:
        rate_ratio = (count_a / exposure_a) / (count_b / exposure_b)
    return {
        "rate_ratio": rate_ratio,
        "exact95_low": convert(lower_probability),
        "exact95_high": convert(upper_probability),
        "p_value": p_value,
    }


def main() -> int:
    args = parse_args()
    with Path(args.opportunities).open(encoding="utf-8") as handle:
        opportunities = list(csv.DictReader(handle, delimiter="\t"))
    with Path(args.candidates).open(encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    eligible_candidates = [row for row in candidates if row.get("parallel_branch_shared", "no") != "yes"]
    excluded_shared = [row["candidate_id"] for row in candidates if row.get("parallel_branch_shared") == "yes"]
    masked_exposure = sum(int(row["site_generations_masked"]) for row in opportunities)
    unmasked_exposure = sum(int(row["site_generations_unmasked"]) for row in opportunities)
    scenarios = [
        (
            "primary_strict_masked_strand",
            "Strict VCF, genotype and BAM validation; mask excluded; bidirectional alternate-read support required.",
            masked_exposure,
            lambda row: row["strict_pre_bam"] == "yes"
            and row["bam_strict_strand_pass"] == "yes"
            and row["in_mask"] == "no",
        ),
        (
            "strict_masked",
            "Strict VCF, genotype and BAM validation; mask excluded; strand requirement relaxed.",
            masked_exposure,
            lambda row: row["strict_pre_bam"] == "yes"
            and row["bam_strict_pass"] == "yes"
            and row["in_mask"] == "no",
        ),
        (
            "strict_unmasked",
            "Strict VCF, genotype and BAM validation without mask or strand exclusion.",
            unmasked_exposure,
            lambda row: row["strict_pre_bam"] == "yes" and row["bam_strict_pass"] == "yes",
        ),
        (
            "relaxed_masked",
            "Relaxed site-quality and allele-balance sensitivity analysis with the mask excluded.",
            masked_exposure,
            lambda row: row["relaxed_pre_bam"] == "yes"
            and row["bam_relaxed_pass"] == "yes"
            and row["in_mask"] == "no",
        ),
    ]
    results: list[dict[str, object]] = []
    selected_ids: dict[str, list[str]] = {}
    selected_by_scenario: dict[str, list[dict[str, str]]] = {}
    for name, description, exposure, keep in scenarios:
        selected = [row for row in eligible_candidates if keep(row)]
        count = len(selected)
        lower_count, upper_count = poisson_interval(count)
        by_interval = Counter(row["interval"] for row in selected)
        results.append(
            {
                "scenario": name,
                "scenario_description": description,
                "mutation_count": count,
                "site_generations": exposure,
                "rate_per_site_per_generation": count / exposure if exposure else 0.0,
                "poisson95_low": lower_count / exposure if exposure else 0.0,
                "poisson95_high": upper_count / exposure if exposure else 0.0,
                "counts_by_interval": json.dumps(dict(by_interval), sort_keys=True),
            }
        )
        selected_ids[name] = [row["candidate_id"] for row in selected]
        selected_by_scenario[name] = selected
    lineage_rates: list[dict[str, object]] = []
    lineage_comparisons: list[dict[str, object]] = []
    lineages = sorted({row.get("lineage", "") for row in opportunities} - {""})
    for name, description, _exposure, _keep in scenarios:
        exposure_field = "site_generations_unmasked" if name == "strict_unmasked" else "site_generations_masked"
        group_values: dict[str, tuple[int, int]] = {}
        for lineage in lineages:
            lineage_intervals = {row["interval"] for row in opportunities if row.get("lineage") == lineage}
            count = len({row["candidate_id"] for row in selected_by_scenario[name] if row["interval"] in lineage_intervals})
            exposure = sum(int(row[exposure_field]) for row in opportunities if row.get("lineage") == lineage)
            lower_count, upper_count = poisson_interval(count)
            group_values[lineage] = (count, exposure)
            lineage_rates.append(
                {
                    "scenario": name,
                    "scenario_description": description,
                    "lineage": lineage,
                    "mutation_count": count,
                    "site_generations": exposure,
                    "rate_per_site_per_generation": count / exposure if exposure else 0.0,
                    "poisson95_low": lower_count / exposure if exposure else 0.0,
                    "poisson95_high": upper_count / exposure if exposure else 0.0,
                }
            )
        if len(lineages) == 2:
            lineage_a, lineage_b = lineages
            count_a, exposure_a = group_values[lineage_a]
            count_b, exposure_b = group_values[lineage_b]
            lineage_comparisons.append(
                {
                    "scenario": name,
                    "lineage_a": lineage_a,
                    "lineage_b": lineage_b,
                    "count_a": count_a,
                    "exposure_a": exposure_a,
                    "count_b": count_b,
                    "exposure_b": exposure_b,
                    **exact_rate_ratio(count_a, exposure_a, count_b, exposure_b),
                }
            )
    output_tsv = Path(args.output_tsv)
    output_json = Path(args.output_json)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    if args.output_lineage_tsv and lineage_rates:
        lineage_path = Path(args.output_lineage_tsv)
        lineage_path.parent.mkdir(parents=True, exist_ok=True)
        with lineage_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(lineage_rates[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(lineage_rates)
    output_json.write_text(
        json.dumps(
            {
                "results": results,
                "selected_candidate_ids": selected_ids,
                "excluded_parallel_branch_shared_candidate_ids": sorted(set(excluded_shared)),
                "lineage_rates": lineage_rates,
                "lineage_comparisons": lineage_comparisons,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
