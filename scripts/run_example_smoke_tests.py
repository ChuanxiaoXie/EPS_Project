#!/usr/bin/env python3
"""Run executable smoke tests against the repository's synthetic fixtures."""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
TESTDATA = REPOSITORY / "testdata"
PYTHON = sys.executable


def run(arguments: list[object], env: dict[str, str] | None = None) -> None:
    command = [str(value) for value in arguments]
    subprocess.run(command, cwd=REPOSITORY, check=True, env=env)


def non_header_records(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def verify_tsv(path: Path, required_columns: set[str], minimum_rows: int = 1) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required_columns.difference(reader.fieldnames or [])
        if missing:
            raise AssertionError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if len(rows) < minimum_rows:
        raise AssertionError(f"{path} contains fewer than {minimum_rows} data rows")


def run_additional_dependency_free_workflows(output: Path) -> dict[str, str]:
    blast_output = output / "filtered_blast.tsv"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "analysis" / "annotation" / "filter_blast_by_coverage.py",
            "--lengths",
            TESTDATA / "annotation" / "query_lengths.tsv",
            "--blast",
            TESTDATA / "annotation" / "blast.outfmt6.tsv",
            "--output",
            blast_output,
        ]
    )
    if len(non_header_records(blast_output)) != 1:
        raise AssertionError("Synthetic BLAST filtering should retain exactly one alignment")

    conservation_output = output / "remapped_conservation.tsv"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "analysis" / "conservation" / "remap_scores.py",
            "--position-map",
            TESTDATA / "conservation" / "position_map.tsv",
            "--scores",
            TESTDATA / "conservation" / "scores.tsv",
            "--output",
            conservation_output,
        ]
    )
    verify_tsv(conservation_output, {"chrom", "start", "end", "score1", "score2"}, 3)

    benchmark_output = output / "variant_benchmark.json"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "analysis" / "simulation" / "benchmark_variants.py",
            "--truth-vcf",
            TESTDATA / "simulation" / "truth.vcf.example",
            "--query-vcf",
            TESTDATA / "simulation" / "query.vcf.example",
            "--output-json",
            benchmark_output,
            "--false-positive-vcf-keys",
            output / "false_positive.tsv",
            "--false-negative-vcf-keys",
            output / "false_negative.tsv",
        ]
    )
    benchmark = json.loads(benchmark_output.read_text(encoding="utf-8"))
    if (benchmark["true_positive"], benchmark["false_positive"], benchmark["false_negative"]) != (2, 1, 1):
        raise AssertionError("Synthetic variant benchmark returned unexpected TP/FP/FN counts")

    junction_output = output / "tdna_junctions.tsv"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "pipelines" / "tdna" / "call_junctions.py",
            "--genome-paf",
            TESTDATA / "tdna" / "genome.paf",
            "--tdna-paf",
            TESTDATA / "tdna" / "tdna.paf",
            "--homology-paf",
            TESTDATA / "tdna" / "homology.paf",
            "--output",
            junction_output,
            "--sample",
            "SyntheticTDNA",
        ]
    )
    verify_tsv(junction_output, {"candidate_id", "contig", "risk_labels"}, 1)

    mutation_status = "skipped (install the mutation-rate optional dependency)"
    if importlib.util.find_spec("scipy") is not None:
        mutation_output = output / "mutation_rate"
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "mutation_rate" / "calculate_rate.py",
                "--opportunities",
                TESTDATA / "mutation_rate" / "opportunities.tsv",
                "--candidates",
                TESTDATA / "mutation_rate" / "candidates.bam_validated.tsv",
                "--output-tsv",
                mutation_output.with_suffix(".tsv"),
                "--output-json",
                mutation_output.with_suffix(".json"),
                "--output-lineage-tsv",
                output / "lineage_rates.tsv",
            ]
        )
        verify_tsv(mutation_output.with_suffix(".tsv"), {"scenario", "mutation_count", "site_generations"}, 4)
        mutation_status = "passed"

    return {
        "annotation_blast_filter": "passed",
        "conservation_remapping": "passed",
        "variant_simulation_benchmark": "passed",
        "tdna_junction_calling": "passed",
        "mutation_rate_calculation": mutation_status,
    }


def render_external_workflows(output: Path) -> None:
    sys.path.insert(0, str(REPOSITORY / "src"))
    from eps_workflows.hic_workflow import render as render_hic
    from eps_workflows.sentieon_gvcf import render_sample
    from eps_workflows.sentieon_joint_calling import (
        read_gvcf_manifest,
        render as render_joint_calling,
    )

    sentieon_config = {
        "reference": str(TESTDATA / "genome" / "reference.fa.example"),
        "sentieon": "/opt/sentieon/bin/sentieon",
        "samtools": "/opt/samtools/bin/samtools",
        "license_server": "license.example.invalid:8990",
        "output_root": "/tmp/synthetic_snp_calling",
        "threads": 2,
        "filter_bam": True,
    }
    sentieon_sample = {
        "sample_id": "SyntheticA",
        "read1": "/data/synthetic_R1.fastq.gz",
        "read2": "/data/synthetic_R2.fastq.gz",
    }
    sentieon_script = render_sample(sentieon_config, sentieon_sample)
    if "Haplotyper" not in sentieon_script or "SyntheticA" not in sentieon_script:
        raise AssertionError("Sentieon renderer did not produce the expected gVCF command")
    if "-F 4 -F 256 -q 20 -f 2 -F 2048" not in sentieon_script or "/^NM:i:/" not in sentieon_script:
        raise AssertionError("Sentieon renderer omitted the post-dedup BAM quality filter")
    (output / "sentieon_gvcf.sh").write_text(sentieon_script, encoding="utf-8")

    joint_script = render_joint_calling(
        {
            "reference": str(TESTDATA / "genome" / "reference.fa.example"),
            "sentieon": "/opt/sentieon/bin/sentieon",
            "gatk": "/opt/gatk/bin/gatk",
            "license_server": "license.example.invalid:8990",
            "output_root": "/tmp/synthetic_joint_calling",
            "joint_vcf_name": "joint.vcf.gz",
            "threads": 2,
        },
        read_gvcf_manifest(TESTDATA / "snp" / "gvcfs.manifest.tsv"),
    )
    if "GVCFtyper" not in joint_script or "VariantFiltration" not in joint_script:
        raise AssertionError("Joint-calling renderer omitted a confirmed publication stage")
    if "haplotype_score > 13.0" not in joint_script:
        raise AssertionError("Joint-calling renderer did not use the reviewed HaplotypeScore direction")
    (output / "sentieon_joint_calling.sh").write_text(joint_script, encoding="utf-8")

    hic_script = render_hic(
        {
            "sample_id": "SyntheticHiC",
            "reference": str(TESTDATA / "genome" / "assembly.fa.example"),
            "read1": "/data/synthetic_hic_R1.fastq.gz",
            "read2": "/data/synthetic_hic_R2.fastq.gz",
            "output_root": "/tmp/synthetic_hic",
            "bwa": "bwa",
            "python": "python3",
            "juicer_root": "/opt/juicer",
            "three_d_dna_root": "/opt/3d-dna",
            "threads": 2,
        }
    )
    if "juicer_v1.sh" not in hic_script or "run-asm-pipeline.sh" not in hic_script:
        raise AssertionError("Hi-C renderer did not produce Juicer and 3D-DNA commands")
    (output / "hic_scaffolding.sh").write_text(hic_script, encoding="utf-8")

    mutation_plan = output / "mutation_rate_plan.json"
    with mutation_plan.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-B",
                str(REPOSITORY / "analysis" / "mutation_rate" / "run.py"),
                "--config",
                str(TESTDATA / "mutation_rate" / "config.json"),
                "--gvcf-manifest",
                str(TESTDATA / "mutation_rate" / "gvcfs.manifest.tsv"),
                "--bam-manifest",
                str(TESTDATA / "mutation_rate" / "bams.manifest.tsv"),
                "--joint-vcf",
                str(TESTDATA / "snp" / "variants.vcf.example"),
                "--reference",
                str(TESTDATA / "genome" / "reference.fa.example"),
                "--mask-bed",
                str(TESTDATA / "mutation_rate" / "mask.bed"),
                "--output",
                str(output / "mutation_rate_pipeline"),
            ],
            cwd=REPOSITORY,
            check=True,
            stdout=handle,
        )
    if len(json.loads(mutation_plan.read_text(encoding="utf-8"))) != 7:
        raise AssertionError("Mutation-rate renderer did not produce the complete command plan")

    tdna_command = subprocess.run(
        [
            PYTHON,
            "-B",
            str(REPOSITORY / "pipelines" / "tdna" / "run_tdna.py"),
            "--config",
            str(TESTDATA / "configs" / "tdna.test.json"),
        ],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if "SyntheticTDNA" not in tdna_command or "pipeline.sh" not in tdna_command:
        raise AssertionError("T-DNA renderer did not produce the expected pipeline command")


def render_top_level_pipelines(output: Path) -> None:
    assembly_config = {
        "assessment_assembly": str(TESTDATA / "genome" / "assembly.fa.example"),
        "output_root": "/tmp/synthetic_assembly_pipeline",
        "python": PYTHON,
        "minimum_sequence_length": 0,
        "merqury": {"enabled": False},
        "busco": {"enabled": False},
    }
    assembly_config_path = output / "assembly_pipeline.json"
    assembly_config_path.write_text(json.dumps(assembly_config, indent=2) + "\n", encoding="utf-8")
    assembly_scripts = output / "assembly_pipeline"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "pipelines" / "genome_assembly" / "run.py",
            "--config",
            assembly_config_path,
            "--script-dir",
            assembly_scripts,
        ]
    )
    assembly_master = assembly_scripts / "run_genome_assembly_pipeline.sh"
    if "assembly_stats.py" not in assembly_master.read_text(encoding="utf-8"):
        raise AssertionError("Genome assembly pipeline omitted assembly statistics")

    annotation_script = output / "genome_annotation.sh"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "pipelines" / "genome_annotation" / "run.py",
            "--config",
            REPOSITORY / "configs" / "annotation.example.json",
            "--script",
            annotation_script,
        ]
    )
    if "repeat_annotation" not in annotation_script.read_text(encoding="utf-8"):
        raise AssertionError("Genome annotation pipeline omitted its first stage")

    gvcf_config = {
        "reference": str(TESTDATA / "genome" / "reference.fa.example"),
        "sentieon": "/opt/sentieon/bin/sentieon",
        "samtools": "/opt/samtools/bin/samtools",
        "license_server": "license.example.invalid:8990",
        "output_root": "/tmp/synthetic_snp_calling",
        "threads": 2,
        "filter_bam": True,
        "samples": [
            {
                "sample_id": "SyntheticA",
                "read1": "/data/synthetic_R1.fastq.gz",
                "read2": "/data/synthetic_R2.fastq.gz",
            }
        ],
    }
    joint_config = {
        "reference": str(TESTDATA / "genome" / "reference.fa.example"),
        "sentieon": "/opt/sentieon/bin/sentieon",
        "gatk": "/opt/gatk/bin/gatk",
        "license_server": "license.example.invalid:8990",
        "output_root": "/tmp/synthetic_joint_calling",
        "joint_vcf_name": "joint.vcf.gz",
        "threads": 2,
    }
    gvcf_config_path = output / "pipeline_gvcf.json"
    joint_config_path = output / "pipeline_joint.json"
    gvcf_config_path.write_text(json.dumps(gvcf_config, indent=2) + "\n", encoding="utf-8")
    joint_config_path.write_text(json.dumps(joint_config, indent=2) + "\n", encoding="utf-8")
    snp_output = output / "snp_pipeline"
    run(
        [
            PYTHON,
            "-B",
            REPOSITORY / "pipelines" / "snp_calling" / "run.py",
            "--gvcf-config",
            gvcf_config_path,
            "--joint-config",
            joint_config_path,
            "--output-dir",
            snp_output,
        ]
    )
    required = [
        snp_output / "samples" / "SyntheticA.sentieon_gvcf.sh",
        snp_output / "generated.gvcfs.tsv",
        snp_output / "sentieon_joint_calling.sh",
        snp_output / "run_snp_calling_pipeline.sh",
    ]
    if any(not path.is_file() for path in required):
        raise AssertionError("SNP-calling pipeline did not render its complete script set")


def locate_rscript() -> str | None:
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
    if os.name == "nt":
        installations = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"), reverse=True)
        if installations:
            return str(installations[0])
    return None


def verify_figure_outputs(paths: list[Path]) -> None:
    for path in paths:
        if not path.is_file() or path.stat().st_size < 1024:
            raise AssertionError(f"R figure output is missing or unexpectedly small: {path}")


def maybe_render_r_figures(output: Path) -> str:
    rscript = locate_rscript()
    if not rscript:
        return "skipped (Rscript not installed)"
    required_packages = "ggplot2,RColorBrewer,circlize,svglite"
    dependency_check = subprocess.run(
        [
            rscript,
            "-e",
            "quit(status=if(all(vapply(c('ggplot2','RColorBrewer','circlize','svglite'), requireNamespace, logical(1), quietly=TRUE))) 0 else 1)",
        ],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
        text=True,
    )
    if dependency_check.returncode != 0:
        return f"skipped (install R packages: {required_packages})"
    fixture_root = output / "r_figure_inputs"
    run([rscript, TESTDATA / "figures" / "generate_r_figure_fixtures.R", fixture_root])
    coverage_fixture = fixture_root / "coverage" / "depth_distribution.tsv"
    coverage_fixture.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(TESTDATA / "coverage" / "depth_distribution.tsv", coverage_fixture)

    hic_prefix = output / "synthetic_hic"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "hic" / "plot_sparse_contact_heatmap.R",
            "--matrix",
            fixture_root / "hic" / "contact_matrix.tsv",
            "--breaks",
            fixture_root / "hic" / "chromosome_breaks.txt",
            "--labels",
            fixture_root / "hic" / "chromosome_labels.txt",
            "--output-pdf",
            hic_prefix.with_suffix(".pdf"),
            "--output-png",
            hic_prefix.with_suffix(".png"),
            "--width",
            "4",
            "--height",
            "4",
            "--dpi",
            "120",
        ]
    )
    snp_prefix = output / "synthetic_snp_circos"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "circos" / "plot_snp_tracks.R",
            "--track-dir", fixture_root / "snp_circos",
            "--chrom-sizes", fixture_root / "snp_circos" / "chrom.sizes.tsv",
            "--output-prefix", snp_prefix,
            "--formats", "pdf,png",
            "--dpi", "120",
        ]
    )
    comparison_prefix = output / "synthetic_comparative_circos"
    comparison_dir = fixture_root / "comparative_circos"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "circos" / "plot_comparative_genome.R",
            "--sectors", comparison_dir / "sectors.tsv",
            "--te", comparison_dir / "te.tsv",
            "--gene", comparison_dir / "gene.tsv",
            "--lai", comparison_dir / "lai.tsv",
            "--links", comparison_dir / "links.tsv",
            "--output-prefix", comparison_prefix,
            "--formats", "pdf,png",
            "--dpi", "120",
        ]
    )
    annotation_prefix = output / "synthetic_annotation_workflow"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "schematics" / "plot_annotation_workflow.R",
            "--output-prefix", annotation_prefix,
            "--formats", "pdf,png",
            "--dpi", "120",
        ]
    )
    snp_workflow_prefix = output / "synthetic_snp_workflow"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "schematics" / "plot_snp_workflow.R",
            "--output-prefix", snp_workflow_prefix,
            "--formats", "pdf,png",
            "--dpi", "120",
        ]
    )
    busco_prefix = output / "synthetic_busco_summary"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "assembly" / "plot_busco_summary.R",
            "--summary", fixture_root / "assembly" / "busco_summary.tsv",
            "--output-prefix", busco_prefix,
            "--formats", "pdf,png",
            "--dpi", "120",
        ]
    )
    coverage_prefix = output / "synthetic_depth_distribution"
    run(
        [
            rscript,
            REPOSITORY / "analysis" / "coverage" / "plot_depth_distribution.R",
            coverage_fixture,
            coverage_prefix.with_suffix(".pdf"),
            coverage_prefix.with_suffix(".png"),
        ]
    )
    prefixes = [
        hic_prefix,
        snp_prefix,
        comparison_prefix,
        annotation_prefix,
        snp_workflow_prefix,
        busco_prefix,
        coverage_prefix,
    ]
    verify_figure_outputs([prefix.with_suffix(extension) for prefix in prefixes for extension in (".pdf", ".png")])
    return "passed (7 R workflows; 14 outputs)"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="eps-public-smoke-") as temporary_directory:
        output = Path(temporary_directory)
        filtered_vcf = output / "hot_region_filtered.vcf"
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "variants" / "hot_region_filter.py",
                "--vcf",
                TESTDATA / "snp" / "variants.vcf.example",
                "--chrom-sizes",
                TESTDATA / "genome" / "chrom.sizes.tsv",
                "--output-vcf",
                filtered_vcf,
                "--hot-regions",
                output / "hot_regions.bed",
                "--window-table",
                output / "window_counts.tsv",
                "--positions-output",
                output / "retained_positions.tsv",
                "--window-size",
                "10",
                "--step-size",
                "5",
                "--hot-count-threshold",
                "2",
            ]
        )
        if len(non_header_records(filtered_vcf)) != 4:
            raise AssertionError("Synthetic hot-region filter should retain exactly four SNPs")

        publication_table = output / "publication_variant_table.tsv"
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "variants" / "build_publication_table.py",
                "--vcf",
                filtered_vcf,
                "--annovar-variant-function",
                TESTDATA / "snp" / "annovar.variant_function.tsv",
                "--annovar-exonic-variant-function",
                TESTDATA / "snp" / "annovar.exonic_variant_function.tsv",
                "--snpeff-tsv",
                TESTDATA / "snp" / "snpeff.tsv",
                "--require-annovar",
                "--output",
                publication_table,
            ]
        )
        verify_tsv(publication_table, {"CHROM", "POS", "Anno_Gene", "SyntheticA_dp_all"}, 4)

        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "variants" / "genotype_table.py",
                "--vcf",
                filtered_vcf,
                "--output-prefix",
                output / "genotypes",
            ]
        )
        verify_tsv(output / "genotypes.sites.tsv", {"chrom", "pos", "heterozygous_samples"}, 4)

        density_input = output / "variant_positions.tsv"
        with density_input.open("w", encoding="utf-8") as handle:
            handle.write("chrom\tpos\n")
            for line in non_header_records(filtered_vcf):
                fields = line.split("\t")
                handle.write(f"{fields[0]}\t{fields[1]}\n")
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "variants" / "window_density.py",
                "--variants",
                density_input,
                "--chrom-sizes",
                TESTDATA / "genome" / "chrom.sizes.tsv",
                "--output",
                output / "window_density.tsv",
                "--window-size",
                "10",
            ]
        )
        verify_tsv(output / "window_density.tsv", {"chrom", "start", "end", "variant_count"}, 1)

        track_dir = output / "circos_tracks"
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "circos" / "prepare_tracks.py",
                "--vcf",
                filtered_vcf,
                "--genes-bed",
                TESTDATA / "genome" / "genes.bed.tsv",
                "--reference",
                TESTDATA / "genome" / "reference.fa.example",
                "--chrom-sizes",
                TESTDATA / "genome" / "chrom.sizes.tsv",
                "--window-size",
                "10",
                "--output-dir",
                track_dir,
            ]
        )
        for filename in ("1_wai.txt", "1_nei.txt", "2_genedensity.txt", "3_GC_wai.txt", "3_GC2AT_nei.txt"):
            if not (track_dir / filename).is_file():
                raise AssertionError(f"Missing synthetic Circos track: {filename}")

        circos_config = json.loads((REPOSITORY / "configs" / "circos.example.json").read_text(encoding="utf-8"))
        circos_config["track_dir"] = track_dir.as_posix()
        circos_config["output_dir"] = (output / "circos_figure").as_posix()
        config_input = output / "circos.test.json"
        config_input.write_text(json.dumps(circos_config, indent=2) + "\n", encoding="utf-8")
        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "circos" / "render_config.py",
                "--config",
                config_input,
                "--output-config",
                output / "circos.conf",
            ]
        )
        if "1_wai.txt" not in (output / "circos.conf").read_text(encoding="utf-8"):
            raise AssertionError("Rendered Circos configuration does not reference the synthetic tracks")

        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "analysis" / "assembly" / "assembly_stats.py",
                "--assembly",
                TESTDATA / "genome" / "assembly.fa.example",
                "--output",
                output / "public_assembly_stats.tsv",
            ]
        )
        verify_tsv(output / "public_assembly_stats.tsv", {"metric", "value"}, 10)

        run(
            [
                PYTHON,
                "-B",
                REPOSITORY / "pipelines" / "tdna" / "fasta_metrics.py",
                TESTDATA / "genome" / "assembly.fa.example",
                output / "assembly_metrics.tsv",
            ]
        )
        verify_tsv(output / "assembly_metrics.tsv", {"metric", "value"}, 4)
        additional_status = run_additional_dependency_free_workflows(output)
        render_external_workflows(output)
        render_top_level_pipelines(output)
        r_status = maybe_render_r_figures(output)

        summary = {
            "dependency_free_workflows": "passed",
            "retained_synthetic_snps": len(non_header_records(filtered_vcf)),
            "circos_tracks": 5,
            "workflow_renderers": ["Sentieon gVCF", "Sentieon joint calling", "Juicer/3D-DNA"],
            "top_level_pipelines": ["genome assembly", "genome annotation", "SNP calling"],
            "r_figure_workflows": r_status,
            "additional_workflows": additional_status,
        }
        print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
