# Input contract v1

The execution plan defines the manifest headers, IDs, paths, counts, and error
policy. This document fixes the previously unspecified analysis settings for
the first validator. Production statistical processing remains P3 work.

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

P1 begins with strict table/sample/count validation; bundle-level design,
reference checksum, and statistical checks must be marked separately until
implemented. The fixture `count_inputs.tsv` is test-only, not a public manifest.

For P0 alignment evidence, the policy matches the legacy gene/exon counting
semantics: unstranded, minimum MAPQ 0, no multimapping inclusion, no overlapping
gene inclusion, no deduplication, PE fragments via `-p --countReadPairs`.
The independent read-origin oracle is SE A=2/B=1/C=0/zero=0, with one ambiguous
read excluded; PE A=1/B=0/C=0/zero=0. These are separate from the hand-authored
technical-run matrix. Both were verified using HISAT2 2.2.1 and Subread 2.0.6.
