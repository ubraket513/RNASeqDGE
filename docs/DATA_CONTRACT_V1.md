# Input contract v1

The execution plan defines the manifest headers, IDs, paths, counts, and error
policy. This document fixes the previously unspecified analysis settings for
the first validator. The P3 offline R command independently enforces the same
contract at its trust boundary.

`analysis.tsv` has unique `key`, `value` rows. Required keys and accepted values:

| Key | Values |
| --- | --- |
| version | `1` |
| design_terms | comma-separated declared column names; `condition` initially |
| alpha | finite decimal strictly between 0 and 1 |
| filter | `zero_total` |
| shrinkage | `apeglm` |

Covariate extensions use `type.<column>` = `categorical` or `numeric`.
`condition` is implicitly categorical. Every other sample column must have an
explicit type. Columns used in a design are identifiers matching
`[A-Za-z][A-Za-z0-9_]*`; terms may not repeat. No arbitrary R expressions are
accepted. Categorical values are nonempty UTF-8 strings, matched exactly.
Numeric values must be finite decimal numbers; missing values are errors.
Numeric predictors cannot be contrast factors. Contrast numerator and
denominator must be distinct observed categorical levels. Their order controls
direction; alphabetical order never chooses a requested comparison.

The reference level for the contrasted factor is its denominator. Other
categorical covariates require `reference.<column>` with an observed level when
they are included in the design. `reference.condition` is optional for a single
condition comparison and, if supplied, must match its denominator. Multiple
contrasts can require separate model releveling in P3. Rank, residual degrees of
freedom, and supported apeglm coefficients require R-side model validation.

### P1 bundle validation decisions

The bundle interface takes five explicit paths, in this order: samples, runs,
references, analysis, contrasts. Each relative data path is resolved against
its own manifest's directory, independently of the current working directory.
The fixture `count_inputs.tsv` is test-only, not a public manifest.

Fixed manifests accept exactly their specified columns (in any order); sample
metadata accepts extra columns only with declared types. Unknown analysis keys,
declarations for absent columns, and duplicate keys are errors. Explicit
`type.condition=categorical` is allowed, but condition cannot be numeric.
Additive designs may use the declared covariates. P3 requires a full-rank model
matrix and positive residual degrees of freedom before fitting.

Every sample must have at least one run. Input paths must resolve to readable
regular files; symlinks are allowed. Paired FASTQ paths must refer to distinct
files. Reusing a file across different run rows is allowed (the synthetic P0
fixture deliberately does this). This validates manifest consistency, not read
content, FASTQ pairing, or biological metadata.

Decimal values accept an optional sign, decimal point and base-10 exponent,
with at least one digit in the significand and exponent when present. Whitespace,
hexadecimal, NaN, infinity and values outside finite double-precision range are
rejected. Alpha additionally lies strictly between zero and one.

References contain exactly the `genome` and `annotation` roles, with nonempty
source and release fields. Digests have exactly 64 hexadecimal digits; letter
case is ignored when comparing SHA-256 against the file contents. These hashes
check the files at validation time; execution snapshots and change detection
remain P5 work.

The denominator supplies the reference for any contrasted categorical factor.
An explicit `reference.<factor>` must agree with all its requested denominators;
omit that fixed reference when multiple contrasts need different denominators.
Every categorical design term not used as a contrast factor requires an explicit
observed reference level. Rank, residual degrees of freedom, releveling and
apeglm coefficient support are validated in P3, not by the native P1 command.
Cross-file count column order and reference gene-set reconciliation belong to P2.

### P2 count adapter and merge decisions

`merge-counts` consumes explicit sample/run manifests, a one-column `gene_id`
reference universe, and an input manifest with `run_id`, `format`, `path`, and
`column`. It supports canonical `gene_id,count` TSV, exact two-column headerless
legacy counts, and featureCounts tables. FeatureCounts selection is explicit:
the first six columns are `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`,
and `column` names the count column. Its adjacent `.summary` file is required,
must contain that same column and a unique `Assigned` row, and all summary/count
values must be unsigned integers. No format or count-column inference occurs.

Every declared run has exactly one count input and every input maps to a declared
run. Count files contain every reference gene exactly once; missing or extra genes
are errors rather than implicit zeros. Inputs may order genes and runs arbitrarily.
Output sorts original gene IDs bytewise and keeps sample columns in `samples.tsv`
order. Technical runs mapped to the same sample are added with checked `uint64`
arithmetic. The aggregate gene-by-sample matrix is retained in memory while each
run file is streamed; P2 therefore claims O(genes × samples + genes) memory, not
constant-memory matrix processing.
The sample ID `gene_id` is reserved because it would collide with the canonical
first output column. Output quoting preserves gene IDs that begin with `#`.

Annotation remains a separate exact-schema table (`gene_id`, `gene_symbol`,
`biotype`, `chromosome`). It may cover a subset of the reference universe and
may contain empty optional values or duplicate symbols, but gene IDs are unique
and must belong to the universe. Annotation validation never joins, drops, or
duplicates count rows.

Original FASTQs may be archived at merge time. The run manifest still validates
IDs, complete sample coverage, layout, strandedness, and required path cells,
while count input existence is checked independently. All merge inputs are
validated before publishing an output. Publication uses a temporary file beside
the destination and atomic rename; failures preserve a previous output and remove
owned temporary files. Output aliases to manifests, count inputs, summaries, or
annotation are rejected. Multi-process locking, resume and provenance snapshots
remain P5 responsibilities.

### P3 offline statistical output decisions

`run_deg_analysis_offline.R` accepts canonical counts, samples, analysis,
contrasts, optional annotation, an absent output path, and a positive worker
count. Counts must be unsigned decimal integers no larger than R's
`.Machine$integer.max`; count columns must exactly equal `samples.tsv` IDs in
their declared order. Only zero-total genes are removed. Annotation is a left
mapping and never changes retained result cardinality.

For every contrast, the denominator is the factor reference. P3 evaluates the
explicit numerator-versus-denominator DESeq2 contrast, verifies it against the
matching named coefficient, and passes that coefficient and raw result to
`lfcShrink(type="apeglm")`. Categorical noncontrast terms use their declared
references; numeric terms must be finite. All input/model checks complete before
the exclusive adjacent stage directory is created.

Each `OUT/contrasts/<contrast_id>/` contains `results.tsv`, `shrunken.tsv`,
`normalized_counts.tsv`, `model_matrix.tsv`, `summary.tsv`, `volcano_data.tsv`,
`pca_data.tsv`, `sample_distances.tsv`, `heatmap_data.tsv`, `sessionInfo.txt`,
and seven PNGs: `MA_raw.png`, `MA_shrunken.png`, `dispersion.png`, `volcano.png`,
`pca.png`, `sample_distances.png`, and `heatmap.png`. Statistical tables retain
every post-filter gene, including NA values. Significance is exactly non-NA
adjusted p-value below `alpha`; PCA and sample distances share the
`varianceStabilizingTransformation(..., blind=FALSE)` assay. The heatmap uses at
most 50 significant genes ordered by adjusted p-value then gene ID. With none,
its TSV has only a header and its PNG explains the empty selection.

The output path must not already exist, including as a symlink, and its parent
must already be writable. Inputs and output may not alias. Successful publication
uses one directory rename; failures remove owned staging data and preserve all
pre-existing paths. Sample IDs `gene_id` and `sample_id`, covariate columns
`PC1`, `PC2`, `PC1_variance`, and `PC2_variance`, and any generated model column
named `sample_id` are reserved to prevent ambiguous output headers.

For P0 alignment evidence, the policy matches the legacy gene/exon counting
semantics: unstranded, minimum MAPQ 0, no multimapping inclusion, no overlapping
gene inclusion, no deduplication, PE fragments via `-p --countReadPairs`.
The independent read-origin oracle is SE A=2/B=1/C=0/zero=0, with one ambiguous
read excluded; PE A=1/B=0/C=0/zero=0. These are separate from the hand-authored
technical-run matrix. Both were verified using HISAT2 2.2.1 and Subread 2.0.6.
