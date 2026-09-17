# Real-data comparison inputs from the report

Sources: `report/main.tex`, `report/Dataset Analysis/dataset_analysis.tex`, and
the GEO SOFT/ENA records pinned in `config/gse80336-sources.json`.
The user confirmed this study is the intended source of real data for P6.

- GEO series: GSE80336 (report section 3.2.4).
- Cohort: 18 bipolar disorder and 18 control samples; also analyze the reported
  35-sample cohort excluding control C_28, GSM2124750. Keep inclusion/exclusion
  explicit in separate manifests; do not infer an outlier removal automatically.
- Covariates described: age, sex, postmortem interval and RIN. The TeX model is
  `~condition` (`main.tex:480,508`); preserve covariates as metadata and use this
  explicit model for the initial report comparison.
- Reference example: GCF_000001405.40, GRCh38.p14, matching FASTA/GTF.
- Run example: SRR3393492. GEO/ENA verify 36 unique experiments/runs, all
  single-end; GEO explicitly describes reverse-stranded libraries.
- The report describes HISAT2 plus HTSeq-count. The new backend comparison uses
  one fixed featureCounts policy for both STAR and HISAT2. Those comparisons must
  be reported separately from reproducing the historical HTSeq-count analysis.
- The report's code requests `control` versus `bipolar`, while nearby prose
  describes bipolar versus control. Declare direction explicitly and reconcile
  it before comparing reported LFC signs. The report's DEG totals also vary
  between sections; do not use one reported number as an unquestioned oracle.

The README previously named BioProject PRJNA316873 alongside GSE80336. Downloaded
GEO SOFT metadata instead explicitly identifies **PRJNA318642 / SRP073382**;
use those verified study relations when staging. The report also mentions GEO supplementary
counts, which may support a real-data statistical check before aligning raw reads.

The GEO SOFT archive was fetched from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80336/soft/GSE80336_family.soft.gz`
and retained under ignored `tests/output/p6-staging/`. It confirms 36 samples,
their titles, age/sex/PMI/RIN metadata and experiment accession relationships.

No raw FASTQ/reference directories were present at the repository's conventional
`data/`, `fastq_files/`, or `ref_genome/` locations when this note was written.
The published count matrix has 47,886 unique Ensembl IDs and 36 samples. It was
processed with TopHat/HTSeq, Ensembl 76, trimming/rRNA removal and upstream gene
symbol collapsing. It is suitable for a statistical pipeline check but cannot
substitute for the STAR/HISAT2 comparison against a common reference/count policy.

`tools/stage_gse80336.R` verifies the downloaded source SHA256 values and writes
explicit full36, without_C28 and balanced subset6 inputs, plus a complete
GSM/SRX/SRR/download-MD5 mapping. All matrices preserve the published integers.
The subset retains three samples per condition, with matched sex distribution
and similar ages; the precise IDs are recorded in the script and sample table.

ENA lists 111,310,546,856 compressed FASTQ bytes. The current machine has about
8 GiB RAM and no `sbatch`. The TeX explicitly identifies Cedar as the historical
compute environment (`main.tex:107,708`). The user now authorizes a reduced local
check on 8 GB RAM, 8 CPUs and 500 GB storage with a 30-minute limit. Stage six
explicit samples with the first 50,000 reads per sample and use chromosome 22 only
(NC_000022.11). This is a restricted-reference workflow check: reads from omitted
chromosomes may map spuriously, and prefix sampling is not random. It does not
establish whole-genome accuracy, performance or biological conclusions. P6 retains
the planned representative subset/full-data, repeated timing and scientific
discrepancy review gates.

## Reproduce the published-count check

Download only the three small study files listed under `files` in
`config/gse80336-sources.json` into `tests/output/p6-staging/`, preserving their
listed filenames. Download is an explicit preparation step; the canonicalizer
and R computation work offline. The reference URLs and NCBI MD5 values are also
pinned in that manifest, but the reference archives have not been downloaded.

```sh
Rscript --vanilla tools/stage_gse80336.R --source tests/output/p6-staging --out tests/output/p6-inputs
Rscript --vanilla tests/integration/check_p6_staging.R
```

For each of `subset6`, `full36`, and `without_C28`, set `cohort` explicitly and
use an absent destination:

```sh
cohort=full36
unshare -Urn env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null \
  .deps/p0/bin/Rscript --vanilla run_deg_analysis_offline.R \
  --counts "tests/output/p6-inputs/$cohort/counts.tsv" \
  --samples "tests/output/p6-inputs/$cohort/samples.tsv" \
  --analysis "tests/output/p6-inputs/$cohort/analysis.tsv" \
  --contrasts "tests/output/p6-inputs/$cohort/contrasts.tsv" \
  --annotation "tests/output/p6-inputs/$cohort/annotation.tsv" \
  --out "tests/output/p6-$cohort" --workers 1
```

The independent oracle is `tests/integration/check_p6_r.R INPUT_DIR OUTPUT_DIR`;
it fits DESeq2/apeglm directly without loading the production analysis code and
compares identities, raw/shrunken statistics, normalization, missing-value masks,
log p-values and DEG membership using the planned tolerances.

## Pending alignment benchmark

Once compute is available, stage FASTQs using the sample/run table URLs and
verify the ENA MD5 values, then verify/decompress the matching NCBI FASTA/GTF.
The current native backend expects uncompressed FASTQ and reference files.
Generate the gene universe from that GTF; do not reuse the published Ensembl
count-matrix IDs as the RefSeq featureCounts universe. Use the same reverse
strand policy, reads, reference, featureCounts settings and `~condition` model
for both backends. Preserve each run's provenance and resource measurements.

Run the explicit six-sample cohort and full 36-sample cohort at least three
times per backend under matched budgets, distinguishing index build/reuse and
cold/warm cache conditions. Record wall/CPU time, peak RSS, disk, mapping and
assignment summaries, count differences, LFC direction/magnitude and DEG sets.
Review discrepancies before any change to the HISAT2 default. The published
counts checks above do not satisfy these pending raw-alignment gates.

## Published-count results

All three cohorts completed with one R worker in a disabled-network namespace:

| Cohort | Input genes | Nonzero genes | padj < 0.05 |
| --- | ---: | ---: | ---: |
| Explicit 3+3 subset | 47,886 | 31,486 | 246 |
| Full 18+18 | 47,886 | 37,540 | 12 |
| Excluding C_28, 18+17 | 47,886 | 37,507 | 15 |

The full-cohort totals agree with the report's 12/15 totals. This alone does not
establish per-gene statistical equivalence. The full/35-sample runs emitted
optimizer warnings; the 35-sample diagnostic log records `nbinomGLM` line-search
failures. Raw and shrunken outputs and warnings are retained separately under
`tests/output/p6-{subset6,full36,without_C28}` and sibling logs. Do not interpret
the much smaller subset's larger DEG count as stronger biological evidence.

The independent six-sample published-count oracle passed raw/shrunken statistics,
normalized counts, gene/sample order, NA/zero masks, log-p values, DEG membership
and seven PNG signatures. An initial checker failure came from apeglm's named
numeric vector attributes; stripping those attributes after separately checking
gene identity resolved it without changing production statistics or tolerances.

## Reduced local run

The 50k-read HISAT2 pilot completed all stages in 92.6 seconds. STAR's first
50k-read sample took 136 seconds to map; the pilot was stopped cleanly to fit
the overall 30-minute limit. Outputs remain under `tests/output/p6-local-benchmark`.
A 10k every-fifth-record attempt under `tests/output/p6-local-10k-benchmark`
failed DESeq2 normalization because every gene had at least one zero. The failure
was retained without a success marker; the statistical method was not changed.
The 50k STAR workflow was resumed, switching its remaining alignment tasks to
two concurrent local workers (two threads each; four requested cores total).
This aims for one functional comparison within the deadline. Three repetitions
and a clean equivalent timing comparison are not claimed. Resource logs include
each attempt, wall/CPU time and maximum process RSS; directory totals exclude
shared indexes. Pauses/restarts and differing task concurrency prevent a clean
end-to-end speed comparison. No OS cache eviction is performed.

Preparation scripts: `tools/stage_gse80336_local.R` fetches/verifies/extracts
chr22 and validates complete FASTQ prefixes; `tools/subsample_gse80336_local.R`
checks prefix hashes before deterministic thinning. ENA HTTPS was unreachable;
the public ENA HTTP endpoint worked. Full reference MD5s are verified, while
partial read downloads cannot verify the full archive MD5 or terminal gzip CRC.
Each local read selection has its own SHA256 and explicit sampling provenance.

`tools/check_gse80336_local.R` runs serial workflows with two native threads,
one R worker, bounded wall time, retained logs and per-run status. The companion
`tools/summarize_gse80336_local.R` requires all requested backend/repetition
pairs before summarizing counts, mapping, junctions, LFC/DEG differences and
repeat consistency. Retain HISAT2 as default regardless of this narrow check.

## Completed local comparison (2026-09-17)

Both 50k-read workflows completed through real DESeq2 and all plots by
13:38:50 UTC, 27 minutes 15 seconds after the 13:11:35 budget start. The independent
summary verifies exact gene/sample order and equality between merged count totals
and featureCounts assigned totals. Evidence is retained in
`tests/output/p6-local-benchmark/comparison.json`; original attempts, local task
commands, per-attempt resource logs, snapshots and scientific outputs are alongside it.

| Metric | HISAT2 | STAR |
| --- | ---: | ---: |
| Input reads (six samples) | 300,000 | 300,000 |
| Mapped primary reads | 8,425 | 60,342 |
| Uniquely mapped primary reads | 7,894 | 51,673 |
| Spliced primary reads | 361 | 3,499 |
| Assigned gene counts | 877 | 985 |
| Nonzero genes | 232 | 267 |
| Significant genes, padj < 0.05 | 0 | 0 |

Across 981 genes and six samples, 105 count cells differ, with summed absolute
difference 126 and count Pearson correlation 0.96209. The largest difference is
RPL23AP82 (17 counts across samples). Among 229 common genes with finite raw LFCs,
correlation is 0.95536, median absolute LFC difference 0.41080, and 14 signs differ.
These sparse-data differences are diagnostics, not biological conclusions.

STAR's much larger mapped-read total yields a smaller assigned-count increase:
most additional alignments do not contribute annotated exon counts. Differences
in mapping/scoring, omitted chromosomes and low-depth sampling prevent an accuracy
ranking. STAR's R fit also reports DESeq2's automatic local dispersion-trend fallback;
the default method was not manually changed. HISAT2 remains the default.

One completed comparison is evidence of functionality, not repeatability. The
three-repeat/full-genome benchmark and live Slurm validation remain deferred under
the user's local constraints. Per-process RSS and interrupted/resumed timings are
retained without presenting them as aggregate peak RAM or a clean speed benchmark.
