#!/usr/bin/env python3
"""Build the manuscript-wide sample genotype table and attach functional annotations."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path
from typing import TextIO

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.variants.genotype_table import classify_genotype, integer_or_none, open_text


SiteKey = tuple[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path, help="Filtered multi-sample VCF.")
    parser.add_argument("--output", required=True, type=Path, help="Wide TSV output; .gz uses gzip.")
    parser.add_argument("--annovar-variant-function", type=Path)
    parser.add_argument("--annovar-exonic-variant-function", type=Path)
    parser.add_argument("--snpeff-tsv", type=Path)
    parser.add_argument(
        "--require-annovar",
        action="store_true",
        help="Drop sites absent from the ANNOVAR variant-function table, matching the legacy join behavior.",
    )
    return parser.parse_args()


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def load_annovar_variant(path: Path | None) -> dict[SiteKey, tuple[str, str]]:
    annotations: dict[SiteKey, tuple[str, str]] = {}
    if path is None:
        return annotations
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: malformed ANNOVAR variant-function row")
            annotations[(fields[2], fields[3])] = (fields[0], fields[1])
    return annotations


def load_annovar_exonic(path: Path | None) -> dict[SiteKey, str]:
    annotations: dict[SiteKey, str] = {}
    if path is None:
        return annotations
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: malformed ANNOVAR exonic-function row")
            protein_change = fields[2].split(":")[-1].rstrip(",") or "NoneChange"
            annotations[(fields[3], fields[4])] = protein_change
    return annotations


def find_column(fieldnames: list[str], candidates: tuple[str, ...]) -> str:
    normalized = {value.lower(): value for value in fieldnames}
    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]
    raise ValueError(f"None of the required columns are present: {', '.join(candidates)}")


def load_snpeff(path: Path | None) -> dict[SiteKey, tuple[str, str]]:
    annotations: dict[SiteKey, tuple[str, str]] = {}
    if path is None:
        return annotations
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        chrom_field = find_column(fieldnames, ("Chr", "CHROM", "chrom"))
        position_field = find_column(fieldnames, ("POS", "pos"))
        effect_field = find_column(fieldnames, ("variant_type", "effect", "annotation"))
        amino_field = find_column(fieldnames, ("Amino_acid_variation", "amino_acid_change", "hgvs_p"))
        for row in reader:
            annotations[(row[chrom_field], row[position_field])] = (row[effect_field], row[amino_field])
    return annotations


def sample_header(sample: str) -> list[str]:
    return [
        sample,
        f"{sample}_dp_all",
        f"{sample}_dp_ref",
        f"{sample}_dp_alt",
        f"{sample}_gq",
        f"{sample}_dp_alt/all",
    ]


def parse_allele_depths(value: str | None) -> tuple[int | None, int | None]:
    if value in {None, "", "."}:
        return None, None
    parsed = [integer_or_none(item) for item in value.split(",")]
    reference_depth = parsed[0] if parsed else None
    alternate_values = parsed[1:]
    alternate_depth = sum(value for value in alternate_values if value is not None) if alternate_values else None
    return reference_depth, alternate_depth


def parse_sample(format_keys: list[str], sample_text: str) -> tuple[list[str], str]:
    call = dict(zip(format_keys, sample_text.split(":")))
    genotype = call.get("GT", "./.").replace("|", "/")
    genotype_class = classify_genotype(genotype)
    if genotype_class == "missing":
        return ["./.", "-", "-", "-", "-", "-"], genotype_class
    depth = integer_or_none(call.get("DP"))
    genotype_quality = integer_or_none(call.get("GQ"))
    reference_depth, alternate_depth = parse_allele_depths(call.get("AD"))
    ratio = alternate_depth / depth if alternate_depth is not None and depth not in {None, 0} else None
    return [
        genotype,
        str(depth) if depth is not None else "-",
        str(reference_depth) if reference_depth is not None else "-",
        str(alternate_depth) if alternate_depth is not None else "-",
        str(genotype_quality) if genotype_quality is not None else "-",
        f"{ratio:.4f}" if ratio is not None else "-",
    ], genotype_class


def build_table(
    vcf_path: Path,
    output_path: Path,
    annovar_variant: dict[SiteKey, tuple[str, str]],
    annovar_exonic: dict[SiteKey, str],
    snpeff: dict[SiteKey, tuple[str, str]],
    require_annovar: bool = False,
) -> dict[str, int]:
    samples: list[str] = []
    written = 0
    dropped_without_annovar = 0
    with open_text(vcf_path) as source, open_output(output_path) as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        for line_number, line in enumerate(source, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                header = [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "Anno_type",
                    "Anno_Gene",
                    "Anno_AAV(Amino acid variation)",
                    "SnpEff_variant_type",
                    "SnpEff_Amino_acid_variation",
                ]
                for sample in samples:
                    header.extend(sample_header(sample))
                header.extend(["ref_sum", "alt_sum", "het_sum", "miss_sum"])
                writer.writerow(header)
                continue
            if not line.strip():
                continue
            if not samples:
                raise ValueError("VCF sample header is missing")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 + len(samples):
                raise ValueError(f"{vcf_path}:{line_number}: too few sample columns")
            chrom, position, _identifier, reference, alternate = fields[:5]
            key = (chrom, position)
            if require_annovar and key not in annovar_variant:
                dropped_without_annovar += 1
                continue
            annotation_type, gene = annovar_variant.get(key, ("NA", "NA"))
            protein_change = annovar_exonic.get(key, "NoneChange" if key in annovar_variant else "NA")
            snpeff_effect, snpeff_amino = snpeff.get(key, ("NA", "NA"))
            row = [
                chrom,
                position,
                reference,
                alternate,
                annotation_type,
                gene,
                protein_change,
                snpeff_effect,
                snpeff_amino,
            ]
            format_keys = fields[8].split(":")
            classes: Counter[str] = Counter()
            for sample_text in fields[9 : 9 + len(samples)]:
                sample_values, genotype_class = parse_sample(format_keys, sample_text)
                row.extend(sample_values)
                classes[genotype_class] += 1
            row.extend(
                [
                    classes["hom_ref"],
                    classes["hom_alt"],
                    classes["heterozygous"],
                    classes["missing"],
                ]
            )
            writer.writerow(row)
            written += 1
    return {"written_sites": written, "dropped_without_annovar": dropped_without_annovar, "sample_count": len(samples)}


def main() -> int:
    args = parse_args()
    counts = build_table(
        args.vcf,
        args.output,
        load_annovar_variant(args.annovar_variant_function),
        load_annovar_exonic(args.annovar_exonic_variant_function),
        load_snpeff(args.snpeff_tsv),
        args.require_annovar,
    )
    print("\t".join(f"{key}={value}" for key, value in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
