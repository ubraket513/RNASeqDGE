# P2 count adapter fixtures

`genes.tsv` lists the four distinct `gene_id` attributes on exon records in
`../p0/reference.gtf`, reviewed by hand. Its reverse order is deliberate: the
merger must emit bytewise sorted IDs. The GTF SHA-256 is
`9bf2a7704737577a810b51df6df0292e2f5abd3c32e399c640dd7128222c0f03`, recorded in
`../p0/references.tsv`. Reference staging must supply this complete gene universe;
P2 does not infer it from a count file or gene symbols.

`inputs.tsv` uses the public P2 schema to point to the existing five P0 normalized
count tables. Their independent expected matrix remains `../p0/expected/counts.tsv`;
it is neither regenerated nor copied from the implementation under test.

`mixed-inputs.tsv` represents the same numbers with shuffled run/gene order and
three formats. The headerless legacy table carries `run_c1a`: A=2, B=3, C=0,
zero=0. The hand-written featureCounts-shaped table carries `run_c1b`: A=4,
B=5, C=1, zero=0. Its unrelated `unused.bam` column contains 99 for every gene
to expose accidental first-column selection. A matching synthetic `.summary`
is retained. These are format/arithmetic fixtures, not generated alignment
results or evidence about alignment quality. The normalized input files are
the existing P0 tables, not independent simulated biological replicates.

The expected `control_1` column is A=2+4=6, B=3+5=8, C=0+1=1, zero=0.
All other sample columns are unchanged. `partial-annotation.tsv` deliberately
omits two genes and duplicates a symbol; supplying it must leave counts identical.

The raw format follows the [Subread users guide, section 6.2.8](https://subread.sourceforge.net/SubreadUsersGuide.pdf)
and was checked against retained local featureCounts 2.0.6 output. Context7 was
consulted first; its result described the R interface, so the official guide
and locally generated pinned-version output supplied the CLI format details.

`make check` compares both adapter paths against the P0 hand-derived matrix.
Optional retained-output validation is a separate explicit command; it uses
`tests/output/p0/reviewed/counts.tsv` and the actual SE/PE featureCounts tables
only when available and never substitutes them for the mandatory oracle.
