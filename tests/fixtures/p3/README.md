# P3 synthetic statistical fixture

Regenerate from the repository root with
`.deps/p0/bin/Rscript --vanilla tools/generate_p3_fixture.R`.
R 4.5.3 seed `17092026` generates 8,000 negative-binomial integers (`size=20`),
with gene means evenly spaced from 40 through 440 and 1,000 rows × 8 columns.
Genes 1–100 are multiplied by five in treated samples; genes 101–200 are
multiplied by five in controls. Batch B columns (2, 4, 6, 8) are multiplied by
two. Genes 201–300 receive an additional twofold treated effect. Genes 501–999
are then reset to zero except for one count in sample `1 + (gene_index %% 8)`.
These sparse rows exercise independent-filtering NA statistics. Gene 1000 is
set to zero in every sample. Each condition has four independent
synthetic replicates, two in each batch. These are synthetic counts, not biological
measurements or evidence of original-study equivalence.

The design is `~ batch + condition`, with batch A as reference. Both contrast
directions are explicit. Annotation covers only genes 1–500 and deliberately
assigns the same symbol to the first two genes. Every original gene ID remains
the row identity. The tests derive their oracle directly with DESeq2, without
sourcing the production script.
