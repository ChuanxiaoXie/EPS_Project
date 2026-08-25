from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from analysis.conservation.remap_scores import parse_old_position
from analysis.circos.prepare_tracks import prepare_tracks
from analysis.circos.render_config import render as render_circos_config
from analysis.mutation_rate.common import bed_bp, load_config, thresholds
from analysis.mutation_rate.calculate_rate import exact_rate_ratio, poisson_interval
from analysis.mutation_rate.identify_candidates import parse_call, target_strict
from analysis.mutation_rate.parse_mpileup import count_bases, pileup_state
from analysis.simulation.benchmark_variants import read_variants, variant_type
from analysis.variants.genotype_table import canonical_substitution, classify_genotype
from analysis.variants.build_publication_table import build_table
from analysis.variants.hot_region_filter import filter_vcf, select_hot_regions
from analysis.variants.window_density import count_variants, read_chromosome_sizes


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


class MutationRateTests(unittest.TestCase):
    def test_public_lineage_configs_are_valid(self) -> None:
        natural = load_config(PACKAGE_ROOT / "configs" / "mutation_rate_natural.example.json")
        transgenic = load_config(PACKAGE_ROOT / "configs" / "mutation_rate_transgenic.example.json")
        self.assertEqual(len(natural["intervals"]), 4)
        self.assertEqual(len(transgenic["intervals"]), 10)
        self.assertEqual(thresholds(natural)["bam_min_mq"], 50)

    def test_call_parsing_and_strict_target_threshold(self) -> None:
        config = load_config(PACKAGE_ROOT / "configs" / "mutation_rate_natural.example.json")
        call = parse_call("0/1:12,8:20:60")
        self.assertEqual(call["code"], 1)
        self.assertAlmostEqual(call["ab"], 0.4)
        self.assertTrue(target_strict(call, thresholds(config)))

    def test_mpileup_base_parser_preserves_strand(self) -> None:
        counts = count_bases(".,Aa", "A")
        self.assertEqual(counts["ref_fwd"], 1)
        self.assertEqual(counts["ref_rev"], 1)
        self.assertEqual(counts["alt_fwd"], 1)
        self.assertEqual(counts["alt_rev"], 1)

    def test_reference_pileup_state(self) -> None:
        config = load_config(PACKAGE_ROOT / "configs" / "mutation_rate_natural.example.json")
        counts = {"ref": 20, "alt": 0, "ab": 0.0}
        self.assertEqual(pileup_state(counts, thresholds(config)), 0)

    def test_bed_opportunity_count(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "regions.bed"
            bed.write_text("Chr1\t0\t10\nChr1\t20\t25\n", encoding="utf-8")
            self.assertEqual(bed_bp(bed), 15)

    @unittest.skipUnless(importlib.util.find_spec("scipy"), "SciPy optional dependency is not installed")
    def test_exact_rate_statistics(self) -> None:
        lower, upper = poisson_interval(0)
        self.assertEqual(lower, 0.0)
        self.assertGreater(upper, 0.0)
        comparison = exact_rate_ratio(2, 100, 1, 100)
        self.assertAlmostEqual(comparison["rate_ratio"], 2.0)
        self.assertGreaterEqual(comparison["p_value"], 0.0)
        self.assertLessEqual(comparison["p_value"], 1.0)


class SupportingAnalysisTests(unittest.TestCase):
    def test_coordinate_identifier_parsing(self) -> None:
        self.assertEqual(parse_old_position("Chr1_123_unused"), ("Chr1", "123"))

    def test_normalized_vcf_key_reader(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            vcf = Path(directory) / "input.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "Chr1\t10\t.\tA\tG\t60\tPASS\t.\n",
                encoding="utf-8",
            )
            variants = read_variants(vcf, {"PASS"})
            self.assertEqual(variants, {("Chr1", 10, "A", "G")})
            self.assertEqual(variant_type(next(iter(variants))), "SNV")

    def test_genotype_and_spectrum_normalization(self) -> None:
        self.assertEqual(classify_genotype("0/1"), "heterozygous")
        self.assertEqual(classify_genotype("1|1"), "hom_alt")
        self.assertEqual(classify_genotype("./."), "missing")
        self.assertEqual(classify_genotype("1/2"), "heterozygous")
        self.assertEqual(canonical_substitution("G", "A"), "C>T")

    def test_hot_region_rule_and_filter(self) -> None:
        windows, regions = select_hot_regions(
            {"Chr1": [9, 19, 29, 149]}, {"Chr1": 200}, 100, 100, 2
        )
        self.assertEqual(regions, {"Chr1": [(0, 100)]})
        self.assertEqual(windows[0], ("Chr1", 0, 100, 3, True))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_vcf = root / "input.vcf"
            output_vcf = root / "output.vcf"
            positions = root / "positions.tsv"
            input_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "Chr1\t10\t.\tA\tG\t60\tPASS\t.\n"
                "Chr1\t20\t.\tC\tT\t60\tLowQual\t.\n"
                "Chr1\t150\t.\tG\tA\t60\tLowQual\t.\n",
                encoding="utf-8",
            )
            counts = filter_vcf(input_vcf, output_vcf, positions, regions, {"PASS"})
            self.assertEqual(counts["retained_snvs"], 2)
            retained = [line for line in output_vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
            self.assertEqual([line.split("\t")[1] for line in retained], ["10", "150"])

    def test_publication_wide_table(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "input.vcf"
            output = root / "table.tsv"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                "Chr1\t10\t.\tC\tT\t60\tPASS\t.\tGT:GQ:DP:AD\t0/1:50:10:6,4\t./.:.:0:0,0\n",
                encoding="utf-8",
            )
            summary = build_table(
                vcf,
                output,
                {("Chr1", "10"): ("exonic", "Gene1")},
                {("Chr1", "10"): "p.A1V"},
                {("Chr1", "10"): ("missense|MODERATE", "p.A1V")},
            )
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(summary["sample_count"], 2)
            self.assertEqual(rows[1][9:15], ["0/1", "10", "6", "4", "50", "0.4000"])
            self.assertEqual(rows[1][-4:], ["0", "0", "1", "1"])

    def test_circos_tracks_and_config(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sizes = root / "sizes.tsv"
            fasta = root / "reference.fa"
            genes = root / "genes.bed"
            vcf = root / "selected.vcf"
            tracks = root / "tracks"
            sizes.write_text("Chr1\t10\n", encoding="utf-8")
            fasta.write_text(">Chr1\nACGTGCAAAT\n", encoding="utf-8")
            genes.write_text("Chr1\t0\t5\tGene1\n", encoding="utf-8")
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "Chr1\t2\t.\tC\tA\t60\tPASS\t.\n"
                "Chr1\t8\t.\tA\tG\t60\tPASS\t.\n",
                encoding="utf-8",
            )
            summary = prepare_tracks(vcf, genes, fasta, sizes, tracks, 5)
            self.assertEqual(summary["gc_to_at_snp_count"], 1)
            with (tracks / "3_GC_wai.txt").open(encoding="utf-8") as handle:
                gc_rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual([row[3] for row in gc_rows], ["3", "1"])
            config_text = render_circos_config(
                {
                    "track_dir": "/data/tracks",
                    "plot_ranges": {
                        "genome_snp": [0, 65],
                        "gene_snp": [0, 25],
                        "gene_density": [0, 15],
                        "gc_count": [0, 119408],
                        "gc_to_at": [0, 49],
                    },
                }
            )
            self.assertIn("file = /data/tracks/1_wai.txt", config_text)
            self.assertIn("max = 119408", config_text)

    def test_variant_density_keeps_groups_separate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sizes_path = root / "sizes.tsv"
            variants_path = root / "variants.tsv"
            sizes_path.write_text("Chr1\t250\n", encoding="utf-8")
            variants_path.write_text(
                "chrom\tpos\tlineage\nChr1\t1\tA\nChr1\t100\tA\nChr1\t101\tB\n",
                encoding="utf-8",
            )
            sizes = read_chromosome_sizes(sizes_path)
            counts, groups = count_variants(
                variants_path, sizes, "chrom", "pos", "lineage", 100
            )
            self.assertEqual(groups, {"A", "B"})
            self.assertEqual(counts[("Chr1", 0, "A")], 2)
            self.assertEqual(counts[("Chr1", 1, "B")], 1)

    def test_all_json_configs_parse(self) -> None:
        for path in (PACKAGE_ROOT / "configs").glob("*.json"):
            with self.subTest(path=path.name):
                self.assertIsInstance(json.loads(path.read_text(encoding="utf-8")), dict)


if __name__ == "__main__":
    unittest.main()
