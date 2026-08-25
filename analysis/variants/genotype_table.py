#!/usr/bin/env python3
"""Convert a multi-sample VCF into auditable site and long-form genotype tables."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--accepted-filter", action="append")
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def integer_or_none(value: str | None) -> int | None:
    if value in {None, "", "."}:
        return None
    return int(float(value))


def classify_genotype(genotype: str) -> str:
    if not genotype or "." in genotype:
        return "missing"
    alleles = genotype.replace("|", "/").split("/")
    if all(value == "0" for value in alleles):
        return "hom_ref"
    if len(set(alleles)) == 1 and alleles[0] != "0":
        return "hom_alt"
    return "heterozygous"


def canonical_substitution(reference: str, alternate: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    if reference in {"A", "G"}:
        reference = reference.translate(complement)
        alternate = alternate.translate(complement)
    return f"{reference}>{alternate}"


def main() -> int:
    args = parse_args()
    accepted_filters = set(args.accepted_filter or ["PASS", "."])
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    site_path = Path(f"{prefix}.sites.tsv")
    genotype_path = Path(f"{prefix}.genotypes.tsv")
    spectrum_path = Path(f"{prefix}.spectrum.tsv")
    site_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "filter",
        "hom_ref_samples",
        "heterozygous_samples",
        "hom_alt_samples",
        "missing_samples",
    ]
    genotype_fields = ["chrom", "pos", "ref", "alt", "sample", "gt", "class", "dp", "ref_ad", "alt_ad", "gq", "alt_fraction"]
    spectrum: Counter[str] = Counter()
    samples: list[str] = []
    with open_text(Path(args.vcf)) as source, site_path.open("w", encoding="utf-8", newline="") as site_handle, genotype_path.open(
        "w", encoding="utf-8", newline=""
    ) as genotype_handle:
        site_writer = csv.DictWriter(site_handle, fieldnames=site_fields, delimiter="\t")
        genotype_writer = csv.DictWriter(genotype_handle, fieldnames=genotype_fields, delimiter="\t")
        site_writer.writeheader()
        genotype_writer.writeheader()
        for line_number, line in enumerate(source, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            if not samples:
                raise ValueError("VCF sample header is missing")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 + len(samples):
                raise ValueError(f"VCF row {line_number} has too few sample columns")
            chromosome, position, _identifier, reference, alternate, _quality, filter_value, _info, format_text = fields[:9]
            if filter_value not in accepted_filters or "," in alternate or alternate in {".", "*"}:
                continue
            format_keys = format_text.split(":")
            classes: Counter[str] = Counter()
            for sample, sample_text in zip(samples, fields[9:]):
                values = sample_text.split(":")
                call = dict(zip(format_keys, values))
                genotype = call.get("GT", ".")
                genotype_class = classify_genotype(genotype)
                classes[genotype_class] += 1
                depth = integer_or_none(call.get("DP"))
                genotype_quality = integer_or_none(call.get("GQ"))
                allele_depths = call.get("AD", ".").split(",")
                reference_depth = integer_or_none(allele_depths[0] if allele_depths else None)
                alternate_depth = integer_or_none(allele_depths[1] if len(allele_depths) > 1 else None)
                denominator = (reference_depth or 0) + (alternate_depth or 0)
                alternate_fraction = alternate_depth / denominator if alternate_depth is not None and denominator else "NA"
                genotype_writer.writerow(
                    {
                        "chrom": chromosome,
                        "pos": position,
                        "ref": reference,
                        "alt": alternate,
                        "sample": sample,
                        "gt": genotype,
                        "class": genotype_class,
                        "dp": depth if depth is not None else "NA",
                        "ref_ad": reference_depth if reference_depth is not None else "NA",
                        "alt_ad": alternate_depth if alternate_depth is not None else "NA",
                        "gq": genotype_quality if genotype_quality is not None else "NA",
                        "alt_fraction": alternate_fraction,
                    }
                )
            site_writer.writerow(
                {
                    "chrom": chromosome,
                    "pos": position,
                    "ref": reference,
                    "alt": alternate,
                    "filter": filter_value,
                    "hom_ref_samples": classes["hom_ref"],
                    "heterozygous_samples": classes["heterozygous"],
                    "hom_alt_samples": classes["hom_alt"],
                    "missing_samples": classes["missing"],
                }
            )
            if len(reference) == len(alternate) == 1 and reference in "ACGT" and alternate in "ACGT":
                spectrum[canonical_substitution(reference, alternate)] += 1
    with spectrum_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["substitution", "site_count"])
        for substitution in ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"):
            writer.writerow([substitution, spectrum[substitution]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
