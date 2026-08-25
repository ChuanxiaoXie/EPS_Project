# Synthetic T-DNA fixtures

The three PAF files form a complete dependency-free test for
`pipelines/tdna/call_junctions.py`. They contain one plant/T-DNA adjacency and
an unrelated homology alignment.

`tdna.fa.example` and `testdata/configs/tdna.test.json` test rendering of the
full assembly pipeline. Executing that full pipeline additionally requires a
real coordinate-sorted BAM plus samtools and the documented assembly tools. A
fake file with a `.bam` suffix is intentionally not supplied because it would
misrepresent integration-test coverage.
