# P0 count-boundary fixture

Hand-authored synthetic data, 2026-09-17. No patient/study metadata, BioMart
responses, aligner output, or DESeq2 estimates are represented here.

Four biological samples have two replicates per condition. The deliberately
nonalphabetical sample order is authoritative. `run_c1a` and `run_c1b` are
technical runs of `control_1`; all other samples have one run. `count_inputs.tsv`
is a **test-only** mapping to normalized count tables, not the production
`runs.tsv` schema or a featureCounts adapter format.

The expected matrix was written independently of any merger. In sample order
treated_2, control_1, treated_1, control_2:

| Gene | Expected counts | Review basis |
| --- | --- | --- |
| gene_A.1 | 20, 6, 18, 8 | control_1 = 2 + 4 |
| gene_B.2 | 4, 8, 2, 9 | control_1 = 3 + 5 |
| gene_C | 7, 1, 6, 2 | control_1 = 0 + 1 |
| gene_zero | 0, 0, 0, 0 | measured zeros in every run |

Gene IDs retain their suffixes. Two IDs share a symbol and another lacks a
symbol; neither case may discard or multiply count rows. Run rows have differing
gene order. All four genes must survive merging; zero-total filtering belongs at
the later statistical boundary. This tiny matrix is not suitable for claiming
DESeq2 dispersion, VST, or plot parity.

`analysis.tsv` instantiates the plan's initial settings: categorical additive
condition design, zero-total filtering, alpha 0.05, apeglm shrinkage. Its exact
keys/values are fixture conventions to carry into P1/P3; a production validator
does not exist yet. The explicit comparison is treated minus control.

Verify the arithmetic and identity contracts from the repository root:

```sh
python3 tests/integration/check_p0_fixture.py
```

This independent stdlib check does not generate or overwrite the expected file.
`expected/` contains a reviewed arithmetic oracle, **not** captured legacy output
or a corrected end-to-end scientific baseline. Legacy observations belong in
`tests/output/p0/legacy.json` and must not be promoted to this oracle.

The sequence fixture includes a deterministic synthetic chromosome with
separate plus/minus genes, a two-exon transcript with a canonical splice junction,
overlapping genes for ambiguous assignment, and a zero-coverage gene; SE reads
and inward-facing PE fragments with explicit origins and strand labels in
`read_origins.tsv`. `references.tsv` records FASTA/GTF SHA-256 hashes; `runs.tsv`
records read paths, layout, and explicit unstranded counting. Recreate these
assets with `python3 tools/generate_p0_sequences.py`. The checker verifies read
sequences against their recorded genomic origins, qualities, and reference hashes.

Runs deliberately reuse reads to test manifests; they are not simulated biological
replicates for statistical validation. The count oracle above is independent of
these reads and **must not be interpreted as expected aligner/counting output**.
The origin labels describe intended cases, not observed mapping or assignment.
Tool-specific tiny-index parameters and assignment expectations still require
review and real HISAT2/STAR/featureCounts execution.
