from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT))
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from analysis.assembly.assembly_stats import SequenceSummary, calculate_metrics
from pipelines.genome_assembly.run import render as render_assembly_pipeline
from eps_workflows.annotation_workflow import render as render_annotation
from eps_workflows.hic_workflow import render as render_hic
from eps_workflows.sentieon_gvcf import render_sample
from eps_workflows.sentieon_joint_calling import render as render_joint_calling


class WorkflowRenderTests(unittest.TestCase):
    def test_assembly_metrics_include_contiguity_and_composition(self) -> None:
        metrics = calculate_metrics(
            [
                SequenceSummary(10, 4, 8, 2),
                SequenceSummary(6, 2, 6, 0),
                SequenceSummary(4, 1, 4, 0),
            ]
        )
        self.assertEqual(metrics["total_bp"], 20)
        self.assertEqual(metrics["N50_bp"], 10)
        self.assertEqual(metrics["L50"], 1)
        self.assertEqual(metrics["N90_bp"], 4)
        self.assertAlmostEqual(metrics["GC_percent_of_ACGT"], 7 / 18 * 100)

    def test_genome_assembly_pipeline_composes_assessment_stages(self) -> None:
        script = render_assembly_pipeline(
            {
                "assessment_assembly": "/data/final.fa",
                "output_root": "/output/assembly",
                "merqury": {
                    "enabled": True,
                    "reads": ["/data/r1.fq.gz", "/data/r2.fq.gz"],
                    "merqury_root": "/opt/merqury",
                },
                "busco": {
                    "enabled": True,
                    "input": "/data/proteins.fa",
                    "lineage": "/db/embryophyta_odb10",
                    "mode": "proteins",
                    "run_name": "BUSCO",
                    "offline": True,
                },
            }
        )
        self.assertIn("assembly_stats.py", script)
        self.assertEqual(script.count("--read"), 2)
        self.assertIn("run_merqury.sh", script)
        self.assertIn("run_busco.sh", script)
        self.assertIn("--mode proteins", script)
        self.assertIn("--offline", script)

    def test_sentieon_renderer_is_parameterized(self) -> None:
        config = {
            "reference": "/reference.fa",
            "sentieon": "/opt/sentieon",
            "license_server": "license:8990",
            "output_root": "/output",
            "threads": 8,
        }
        sample = {"sample_id": "SAMPLE01", "read1": "/r1.fq.gz", "read2": "/r2.fq.gz"}
        script = render_sample(config, sample)
        self.assertIn("SAMPLE01.g.vcf.gz", script)
        self.assertIn("--algo Haplotyper", script)
        self.assertNotIn("rm -f", script)

    def test_sentieon_renderer_includes_optional_post_dedup_filter(self) -> None:
        config = {
            "reference": "/reference.fa",
            "sentieon": "/opt/sentieon",
            "samtools": "/opt/samtools",
            "license_server": "license:8990",
            "output_root": "/output",
            "threads": 8,
            "filter_bam": True,
            "minimum_mapping_quality": 20,
            "maximum_edit_distance": 1,
            "maximum_cigar_operations": 2,
        }
        sample = {"sample_id": "SAMPLE01", "read1": "/r1.fq.gz", "read2": "/r2.fq.gz"}
        script = render_sample(config, sample)
        self.assertIn("-F 4 -F 256 -q 20 -f 2 -F 2048", script)
        self.assertIn("/^NM:i:/", script)
        self.assertIn("max_nm=1", script)
        self.assertIn("max_ops=2", script)
        self.assertIn("SAMPLE01.rmdup.filtered.bam", script)

    def test_hic_renderer_is_non_destructive(self) -> None:
        config = {
            "sample_id": "EPS",
            "reference": "/reference.fa",
            "read1": "/r1.fq.gz",
            "read2": "/r2.fq.gz",
            "output_root": "/hic",
            "bwa": "bwa",
            "python": "python3",
            "juicer_root": "/opt/juicer",
            "three_d_dna_root": "/opt/3d-dna",
        }
        script = render_hic(config)
        self.assertIn("safe_link", script)
        self.assertIn("merged_nodups.txt", script)
        self.assertNotIn("rm -rf", script)

    def test_joint_calling_renderer_reproduces_confirmed_stages(self) -> None:
        config = {
            "reference": "/reference.fa",
            "sentieon": "/opt/sentieon",
            "gatk": "/opt/gatk",
            "license_server": "license:8990",
            "output_root": "/output",
            "joint_vcf_name": "joint.vcf.gz",
        }
        script = render_joint_calling(config, ["/input/A.g.vcf.gz", "/input/B.g.vcf.gz"])
        self.assertIn("--algo GVCFtyper", script)
        self.assertIn("SelectVariants", script)
        self.assertIn("VariantFiltration", script)
        self.assertIn("QD < 2.0", script)
        self.assertIn("haplotype_score > 13.0", script)
        self.assertEqual(script.count("-v /input/"), 2)

    def test_annotation_stages_are_dependency_ordered(self) -> None:
        config = {
            "work_dir": "/annotation",
            "stages": [
                {"name": "integrate", "depends_on": ["align"], "command": "run_integrate"},
                {"name": "align", "command": "run_align"},
            ],
        }
        script = render_annotation(config)
        self.assertLess(script.index("run_stage align"), script.index("run_stage integrate"))

    def test_tdna_pipeline_is_parameterized(self) -> None:
        pipeline = (PACKAGE_ROOT / "pipelines" / "tdna" / "pipeline.sh").read_text(encoding="utf-8")
        self.assertIn("Set BAM_PRIMARY", pipeline)
        self.assertIn("SCRIPT_DIR", pipeline)

    def test_example_configs_are_valid_json(self) -> None:
        for path in (PACKAGE_ROOT / "configs").glob("*.json"):
            with self.subTest(path=path.name):
                with path.open(encoding="utf-8") as handle:
                    self.assertIsInstance(json.load(handle), dict)


if __name__ == "__main__":
    unittest.main()
