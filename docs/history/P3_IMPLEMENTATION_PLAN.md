# P3 offline DESeq2 implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development to implement this plan with test-first
> changes and independent review.

**Goal:** Add a pinned, network-free R/DESeq2 analysis interface over the P2
canonical counts and v1 manifests, with explicit designs, contrasts, references,
auditable plot data, and numerical parity tests.

**Architecture:** Add `run_deg_analysis_offline.R` as the P3 entry point. Keep
the existing GEO-driven `run_deg_analysis.R` unchanged until the later workflow
integration gate, so the active legacy pipeline retains its recovery path. The
new script validates inputs itself, builds additive model formulas only from
declared identifiers, fits once per explicit contrast/reference, and writes each
contrast into an isolated output directory.

**Tech stack:** R 4.5.3, DESeq2 1.50.2, apeglm 1.32.0, optparse,
BiocParallel, ggplot2, pheatmap, RColorBrewer. No network, GEOquery, BioMart,
org.Hs.eg.db, shell commands, or runtime package installation.

**Spec:** `docs/history/IMPLEMENTATION_PLAN.md` P3/statistical boundary and
`docs/DATA_CONTRACT_V1.md`.

## Global constraints

- Inputs are canonical TSV: counts, samples, analysis, contrasts, and optional
  annotation. Counts/sample identity and order are exact; no intersection.
- Accept only nonnegative decimal integer counts no larger than
  `.Machine$integer.max`; reject silent rounding/coercion.
- `filter=zero_total`, `shrinkage=apeglm`, and explicit `alpha` are binding.
- Design is additive declared terms only. Categorical/numeric typing and
  references follow v1. Check finite numeric predictors, full model rank, and
  positive residual degrees of freedom before fitting.
- For each contrast, set the contrasted factor reference to its denominator,
  use `results(..., contrast=c(factor,numerator,denominator), alpha=alpha)`,
  identify and verify the matching coefficient, then pass that raw result and
  coefficient to `lfcShrink(type="apeglm")`. Direction must never be inferred
  alphabetically.
- Each contrast writes `results.tsv`, `shrunken.tsv`, `normalized_counts.tsv`,
  `model_matrix.tsv`, `summary.tsv`, `volcano_data.tsv`, `pca_data.tsv`,
  `sample_distances.tsv`, `heatmap_data.tsv`, seven nonempty PNGs, and
  `sessionInfo.txt` under `OUT/contrasts/<contrast_id>/`.
- Result tables retain every post-filter gene, including NA statistics, and use
  `gene_id` as the first column. Optional annotation is a left mapping only;
  missing symbols stay missing and duplicate symbols never change row count.
- Significant membership is exactly non-NA `padj < alpha` in results, summaries,
  volcano data, and significant heatmaps. Volcano x is shrunken LFC from the same
  contrast; y/significance use its verified raw result. Zero significant genes
  is valid: write header-only heatmap data and an explanatory nonempty PNG.
- Use `varianceStabilizingTransformation(..., blind=FALSE)`. PCA and sample
  distances use the same transformed assay; heatmap rows use significant genes
  ordered by padj, capped at 50. Handle one-row heatmaps without invalid
  clustering.
- `--workers 1` uses `SerialParam`; larger values use `MulticoreParam`. Outputs
  from workers 1/2/4 must have identical row/sample membership and numerical
  agreement `abs(a-b) <= 1e-10 + 1e-7*abs(reference)`, with matching NA masks.
- Validate all inputs and models before creating result files. Write into an
  exclusive temporary output directory beside `OUT`, then atomically rename;
  preserve an existing `OUT` on any failure and reject input/output aliasing.
- P3 claims parity only against an independent direct-DESeq2 oracle on a synthetic
  fixture in the same pinned environment. It does not claim original-study parity.

## Ruling

Add a new offline script rather than replacing `run_deg_analysis.R`. The active
Snakemake workflow still supplies `--gse`; replacing the file during P3 would
break the preserved HISAT2 baseline before P5 integration. Cost if wrong: P5 must
explicitly switch to the offline command instead of inheriting the filename.

## Task 1: offline backend and scientific fixture

**Files:**
- Create: `run_deg_analysis_offline.R`
- Create: `tools/generate_p3_fixture.R`
- Create: `tests/fixtures/p3/{counts,samples,analysis,contrasts,annotation}.tsv`
- Create: `tests/integration/check_p3_r.R`

**Interfaces:**
- Command: `Rscript --vanilla run_deg_analysis_offline.R --counts FILE
  --samples FILE --analysis FILE --contrasts FILE [--annotation FILE]
  --out DIR [--workers N]`.
- Output paths and semantics are exactly those in Global constraints.

- [ ] Generate a deterministic 1,000-gene/eight-sample fixture with balanced
  `batch` and `condition`, at least three biological replicates per condition,
  known treated-up genes, one zero-total gene, a partial annotation with duplicate
  symbols, and forward/reverse condition contrasts. Record the seed and arithmetic
  construction in the fixture README.
- [ ] Write integration tests first. The positive test invokes the absent command
  and must fail for the missing offline script. Add independent direct-DESeq2
  calculations for the forward contrast and literal checks for input/sample/gene
  identities, direction, NA masks, result tolerance, annotation row preservation,
  contrast inversion, plot-data semantics, PNG nonemptiness, and session info.
- [ ] Implement strict manifest/count/model validation and the per-contrast fit,
  output tables, plots, and atomic directory publication.
- [ ] Add negative cases for reordered/missing samples, count overflow/fraction,
  undeclared/missing covariates, numeric NA/nonfinite, absent contrast levels,
  rank deficiency, zero residual df, unknown annotation genes, output aliases,
  and preservation of a pre-existing output directory.
- [ ] Run workers 1, 2, and 4; compare membership, NA masks, and numeric tolerance.
  Run a separate no-signal fixture and verify the zero-significant output path.
- [ ] Run the pinned R integration test and the existing `make -j2 check` gate.
  Write a concise report; do not commit because P1/P2 are intentionally
  uncommitted in the user-authorized shared `main` checkout.

## Task 2: review, documentation, and completion gate

**Files:**
- Modify: `docs/DATA_CONTRACT_V1.md`
- Modify: `docs/LOCAL_DEVELOPMENT.md`
- Modify: `docs/IMPLEMENTATION_PROGRESS.md`
- Modify: `docs/history/SESSION_HANDOFF.md`
- Modify: Serena project memories if stable continuation facts changed.

- [ ] Independently review the R interface for statistical direction, coefficient
  matching, rank/DF checks, sample identity, NA handling, FDR consistency, output
  publication, test independence, and legacy-path preservation.
- [ ] Fix all material findings with regression tests and scoped re-review.
- [ ] Run fresh `make -j2 check`, pinned `check_p3_r.R`, workers 1/2/4 parity,
  `git diff --check`, and a source scan confirming the offline script has no GEO,
  BioMart, download, install, or shell execution.
- [ ] Record exact results and remaining P4/P5/P6 limits. Do not report runtime
  estimates, original-study parity, or backend equivalence as measured evidence.

## Preflight interface review

| Tasks | Shared surface | Finding / resolution |
| --- | --- | --- |
| 1 / legacy workflow | R script filename and CLI | use a new offline file; P5 owns workflow switch |
| 1 / P2 | canonical count/sample/annotation schemas | exact reuse; no sample intersection or annotation join filtering |
| 1 / P1 | analysis/contrast validation | repeat at R trust boundary; P1 remains early native validation |
| 1 / 2 | output schema and parity evidence | fixed in Global constraints; docs report only executed gates |
| Task 1 internal | multiple denominators vs apeglm coefficient | fit each contrast after denominator releveling; verify raw contrast equals matched coefficient |
| Task 1 internal | significant heatmap vs zero DEG | significant-only, max 50; explicit empty-data/PNG behavior |
