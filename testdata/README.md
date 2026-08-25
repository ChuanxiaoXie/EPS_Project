# Synthetic smoke-test data

This directory contains small, fully synthetic inputs for executable examples.
They contain no biological observations, real sample identifiers or private
storage paths and must not be used for scientific interpretation.

The fixtures cover:

- assembly FASTA metrics and workflow rendering;
- Hi-C sparse contact-heatmap inputs;
- multi-sample SNP filtering and table generation;
- Sentieon joint-genotyping and GATK hard-filter script rendering;
- five-track Circos data preparation; and
- ANNOVAR- and SnpEff-shaped annotation joins;
- BLAST filtering, conservation remapping and simulation benchmarking;
- coverage plotting, mutation-rate calculation and T-DNA junction calling;
- GO-enrichment input shapes; and
- protein-structure input for the optional SaProt integration.

`testdata/test_matrix.tsv` maps every public analysis, pipeline or workflow
entry point to its fixtures and test mode. Run
`python scripts/check_testdata_coverage.py` to verify that every declared
fixture exists and that no public entry point is omitted.

Test modes have precise meanings: `executed` means the analysis runs against
synthetic data; `rendered` means its command or shell script is generated and
inspected; `unit_tested` covers internal parser logic; and
`optional_integration` identifies an external licensed tool, model, database or
specialist executable needed for a complete run.

`testdata/figures/generate_r_figure_fixtures.R` creates deterministic
ten-chromosome inputs for the Hi-C, SNP Circos and comparative-genome R plots.
The generated files are temporary smoke-test artifacts and are not manuscript
source data.

`testdata/assembly/busco_summary.tsv` supplies artificial BUSCO counts for the
assembly-summary plotting test.

Run all dependency-free smoke tests from the repository root:

```bash
python scripts/run_example_smoke_tests.py
```

If `Rscript` with ggplot2, RColorBrewer, circlize and svglite is available, the
same command renders PDF and PNG outputs for all seven core R figure workflows.
Sentieon, Juicer, 3D-DNA, Merqury, Meryl, BUSCO lineage databases, SaProt model
weights and Foldseek are not bundled; their commands or interfaces are checked
without pretending those external programs are installed. GO enrichment also
requires the Bioconductor packages `goseq` and `topGO`.

All fixtures are artificial and deliberately tiny. They test file formats,
parameter flow and expected control paths, not biological validity.
