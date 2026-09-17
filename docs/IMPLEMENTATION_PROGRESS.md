# Implementation progress

Updated 2026-09-17. Latest user checkpoint: `2569861` (performance profiling).

## Status

**P0 local gate passed; P1 validation, P2 streaming count merge, and P3 offline
DESeq2 interface and P4 offline alignment/toolchain boundary complete.
P5 complete; P6 reduced local comparison complete; full benchmark deferred;
P7 complete for the authorized local scope.** Native scheduling/profiling and
compatible R fit reuse passed their checks. The Python-free compute prefixes
passed package/library audits, eight real-tool cases and a six-sample raw-read
workflow with exact merged counts and all nine scientific tables matching P6.
See [P7 runtime evidence](P7_RUNTIME.md). Helper retirement, sampling-provenance
review fix and Slurm submission-affinity regression are complete. HISAT2 remains the default. Dependencies
were installed after explicit user authorization. Production entry points now
use the native workflow; no live Slurm jobs were submitted. Published-count checks and a
reduced real-read workflow have run; whole-study raw-alignment parity is not claimed.
The pre-existing README changes and untracked docs/Serena metadata were preserved.

## Continuation: dependency installation and implementation

The initial-session sections below are retained as historical evidence. The
missing-dependency blockers there are resolved for the local checks described
here. SRA Toolkit and Slurm are still absent; no remote downloads/jobs are needed
for these tests.

### Installed environments

Used Mamba 2.5.0 with conda-forge then bioconda and strict channel priority.
The initial combined solve failed because samtools 1.17's old zlib requirement
conflicted with R. Installed isolated project prefixes instead:

- `.deps/p0`: 344 packages, Python 3.11, Snakemake 8.30.0, R 4.5.3,
  DESeq2 1.50.2, apeglm 1.32.0, and all 12 required R libraries.
- `.deps/alignment`: 48 packages, HISAT2 2.2.1, STAR 2.7.10b,
  Subread/featureCounts 2.0.6, samtools 1.17.

Explicit package URL/build/MD5 locks and SHA-256/license/dependency provenance
are in `config/*-linux-64.*`. All 392 conda package records have SHA-256 hashes.
The org.Hs.eg.db post-link installer separately downloaded its source archive,
verified MD5, and removed it. Its URLs/MD5 are recorded from the pinned installed
recipe; an extra source-archive SHA-256 was not captured. These are online
staging environments, not an offline redistribution bundle or P4 source-build
manifest. See [LOCAL_DEVELOPMENT.md](LOCAL_DEVELOPMENT.md) for reproduction.

### P0 execution evidence

| Check | Actual result |
| --- | --- |
| `python tests/integration/check_p0_fixture.py` | 3 tests pass |
| `python tests/integration/check_p0_legacy.py` in installed environment | 4 tests pass; reproduces gene loss, missing-sample skip, and unaggregated technical runs; corrected gene-ID reference matches independent oracle |
| `python tests/integration/check_p0_alignment.py` | SE and PE pass with real HISAT2/featureCounts; BAM integrity, exact gene counts and assignment summaries checked; junction CIGAR `37M200N38M` verified |
| `Rscript --vanilla tests/integration/check_p0_r.R` | All libraries load; 1,000-gene/6-sample synthetic DESeq2 and apeglm calculation, explicit coefficient/contrast, and PNG pass; sessionInfo saved |
| Original Snakemake default dry-run | Only `make_directories` scheduled, confirming default-target issue |
| Original awk/R-expression probe after install | Merger import now works; tab and suffix-removal probes still pass |

Original merger tables are under `tests/output/p0/original/`; the reviewed
gene-ID merged matrix is under `tests/output/p0/reviewed/`. Neither is conflated
with original-study output. The hand-authored matrix remains unchanged.
All alignment invocations and raw featureCounts/summary files are retained under
`tests/output/p0/alignment/`. Local tools run sequentially with one worker;
samtools sort uses zero extra threads. No multithread invariance claim is made.

R logs record a parametric-to-local dispersion fit fallback on the synthetic
data. This is an environment functionality test, not statistical migration
parity. The harmless timedatectl warning reflects the sandbox's unavailable
system bus; R completed successfully. Full dataset/cluster validation remains
P6 work. [DATA_CONTRACT_V1.md](DATA_CONTRACT_V1.md) fixes analysis keys and the
proposed covariate declarations and documents the independently derived read
assignment oracle.

### P1 implementation and verification

- Added C++20 `src/`, `include/rnaseq/`, and Make build with user compiler/linker
  flags and generated dependencies, including vendored headers. No download,
  R, Python, Slurm, host-wide worker detection, or parser threads in native checks.
- Vendored csv-parser 5.3.0 release header, verified against upstream SHA-256,
  and doctest 2.4.12 from a pinned commit. Licenses, embedded dependency notices,
  and checksum manifest are in `third_party/`. Context7 returned unrelated
  libraries, so the exact upstream release source/header was inspected.
- Added strict UTF-8/TSV syntax validation, sample/count schema commands,
  duplicate checks, checked uint64 parsing/addition, and location-bearing errors.
  Quoted fields, BOM and CRLF are covered; malformed fields and missing files
  fail. Output stream failures are tested against `/dev/full`.
- `make -j2 check`: **6 test cases / 58 assertions pass**, CLI checks pass.
- `make -j2 -B sanitize`: AddressSanitizer/UndefinedBehaviorSanitizer and leak
  checks pass, including CLI checks. The first sandbox run hit LeakSanitizer's
  ptrace restriction; the approved rerun outside the sandbox succeeded.
- `make -n -W third_party/csv-parser/csv.hpp all` confirms vendor header changes
  schedule recompilation. GCC 15.2.0 was tested; GCC 12 is still a target minimum,
  not a tested compatibility claim.

At the end of that initial implementation, P1 was not complete: `validate-tsv`
was syntax-only for run/reference/analysis/contrast manifests. The continuation
below completes their semantic/cross-file validator. The current
table loader materializes inputs in memory and is not a large-count performance
implementation. R integer limits/model validation and statistical parity remain
P3. CPU STAR is installed but has not yet been validated as an alignment backend.

### P1 continuation (2026-09-17)

Resumed from committed `f2c33ca` with a clean worktree. User explicitly requested
working in the existing checkout on `main`. Initial `make -j2 check` passed the
existing 6 test cases / 58 assertions and CLI checks. Serena memories and text
search work; C++ symbol navigation is unavailable because the project is
configured for Bash only.

Completed: native bundle semantics, safe reference hashing, unit/CLI regression
coverage, combined native and sanitizer verification, and independent review. The
manifest implementation and reference hashing have separate source/test files;
their shared interface is `sha256_file(path)`, which returns a digest or throws.
The CLI and Makefile integrate both. No P2/P3 behavior is included in this gate.

Contract decisions are recorded in `DATA_CONTRACT_V1.md`: strict fixed schemas
and known settings; complete sample/run coverage; finite decimal/exponent
syntax without whitespace; readable regular files with symlinks allowed;
distinct paired files; case-insensitive SHA-256; explicit covariate references
consistent with requested denominators. These resolve previously unspecified
input cases and may reject manifests that only passed syntax validation.

Reference hashing reuses the existing `sha256sum` dependency through POSIX
`posix_spawnp`, with an opened input descriptor and bounded output, without a
shell. Context7 GNU coreutils documentation confirmed standard-input handling.
Local `sha256sum --version` reports uutils coreutils 0.8.0. Known empty/abc/
million-a vectors, unusual filenames, missing/nonregular files, and utility
failures/malformed output are covered. The initial vector tests failed with the
unimplemented stub; the final focused hashing suite passed 4 cases / 31 assertions.

Added `include/rnaseq/{manifests,file_hash}.hpp`, `src/{manifests,file_hash}.cpp`,
and their unit tests; integrated them into `src/main.cpp`, `Makefile`, and the
native CLI checks. `validate-bundle SAMPLES RUNS REFERENCES ANALYSIS CONTRASTS`
validates the complete P0 fixture. Existing single-table commands retain their
scope. No source data is rewritten, external analysis tool run, or job submitted.

Review found and fixed bare-CLI argument handling, symlink/`..` path resolution,
and an overly strict identifier restriction on metadata-only columns. Regression
tests cover these fixes, absolute paths, paired-file equivalence, multiple
contrasts, numeric syntax, missing mappings/files, mismatched reference hashes,
and invalid design/contrast settings. The focused re-review found no new issues.
Future multithreaded subprocess execution and general timeout/signal policy
remain P5 work; P1 uses a single thread and a trusted checksum utility on PATH.

Final verification on the combined implementation:

| Command | Actual result |
| --- | --- |
| `make -j2 check` | 18 test cases / 131 assertions pass; vendor checks and CLI integration pass |
| `make -j2 sanitize` | ASan/UBSan and leak checks pass; all 18 cases / 131 assertions and CLI integration pass |
| `git diff --check` | Pass |

The sandbox sanitizer attempt passed assertions but hit LeakSanitizer's ptrace
restriction; the explicitly approved rerun outside the sandbox passed in full.
No new dependency installation or statistical-parity claim was made. Changes
remain uncommitted in the existing checkout as requested. Updated Serena's stale
core/completion memories to point to these documents and the native test suite.

### P2 streaming count adapters and merge (2026-09-17)

Added a reusable `TsvReader` that validates and decodes one physical record at a
time with the pinned csv-parser, without normalizing the complete input in memory.
Small manifests remain materialized. Added canonical TSV, exact two-column legacy,
and featureCounts adapters with explicit column selection and required `.summary`
validation. All declared run inputs contain the complete explicit reference gene
universe. Technical runs aggregate with checked `uint64` arithmetic into a matrix
ordered by bytewise gene ID and authoritative sample order.

Added `merge-counts`, `validate-counts-for-samples`, and `validate-annotation`.
Annotation remains an independently validated subset and cannot filter or multiply
counts. Merge permits archived FASTQs while preserving run/sample/layout/strand
validation. Inputs are validated before an exclusive adjacent temporary output is
atomically renamed; failures preserve prior output and aliases to input files are
rejected. Multi-process locking, resume, and provenance snapshots remain P5.

The reproducible P2 fixtures cover all three formats, shuffled genes/runs, a
distractor featureCounts column, technical-run aggregation, partial annotation,
overflow, missing/extra genes, summaries, path semantics, nonregular files, and
failure-safe publication. Both canonical and mixed-format merges match the
unchanged hand-derived P0 matrix byte-for-byte. The optional retained-evidence
check also matches `tests/output/p0/reviewed/counts.tsv` and the independently
reviewed real featureCounts 2.0.6 SE/PE assignment oracle; this is not an
alignment-backend or statistical-parity claim.

Independent review found two output-boundary defects: leading-`#` gene IDs were
not re-quoted, and sample ID `gene_id` collided with the canonical first header.
Both now have regression tests; the scoped re-review found no new breakage.

| Command | Actual result |
| --- | --- |
| `make -j2 check` | 34 test cases / 215 assertions; native and count CLI integrations pass |
| `bash tests/integration/check_p2_retained.sh build/rnaseq` | reviewed matrix and retained real SE/PE featureCounts outputs pass |
| `make -j2 sanitize` | ASan/UBSan/leak checks pass with 34 cases / 215 assertions and both CLI integrations |
| `git diff --check` | Pass |

No production pipeline, backend default, R analysis, scheduler state, or original
study output changed. The current aggregate uses O(genes × samples + genes) memory;
P2 does not claim constant-memory output, concurrency safety, or performance gains.

### P3 offline R/DESeq2 interface (2026-09-17)

Added `run_deg_analysis_offline.R`, a network-free interface over the canonical
counts, samples, analysis, contrasts, and optional annotation TSVs. It repeats
strict validation at the R boundary, preserves sample and retained gene identity,
checks full model rank and positive residual degrees of freedom, and fits each
requested direction after releveling its stated denominator. Raw DESeq2 contrast
results are verified against the coefficient used by apeglm before shrinkage.

Each contrast is staged beside the requested absent output directory and then
published by one rename. Existing outputs and input aliases are rejected. Results
retain NA statistics and use one `padj < alpha` membership rule for summaries,
volcano data, and significant heatmaps. A zero-significant result is valid and
produces header-only heatmap data plus an explanatory PNG. Reserved sample and
covariate names that would duplicate generated TSV headers are rejected.

The deterministic 1,000-gene/eight-sample fixture is checked through the public
CLI against separately constructed DESeq2/apeglm calculations. It covers forward
and reverse contrasts, workers 1/2/4, NA masks, annotation cardinality, numeric
covariates, absent annotation, unusual quoted IDs, no- and one-significant-gene
heatmaps, malformed inputs/models, fitting failure, output preservation, and all
seven PNGs per contrast. The legacy `run_deg_analysis.R` and active Snakemake
workflow remain unchanged for the later P5 switch.

P3 establishes only synthetic same-environment DESeq2 parity. It does not
establish original-study parity, alternative-backend equivalence, workflow
integration, caching/resume, provenance snapshots, concurrent publication
locking, or performance improvement.

Fresh completion verification after the reserved-output-column regression fix:

| Command | Actual result |
| --- | --- |
| `.deps/p0/bin/Rscript --vanilla tests/integration/check_p3_r.R` | Pass; independent forward DESeq2/apeglm and reverse DESeq2 oracles, workers 1/2/4 parity, 26 malformed/model/fitting cases, alias and existing-output preservation, all artifact/heatmap paths |
| `make -j2 check` | 34 test cases / 215 assertions; native and count CLI integrations pass |
| `git diff --check` | Pass |
| offline source scan | No GEOquery, BioMart, download, install, shell/system, source, or URL calls in `run_deg_analysis_offline.R` |
| fixture regeneration/hash comparison | All five P3 TSV hashes unchanged |

The pinned analysis versions observed by the suite are R 4.5.3, DESeq2 1.50.2,
apeglm 1.32.0, BiocParallel 1.44.0, and pheatmap 1.0.13. Logs are retained under
ignored `tests/output/p3/`.

### P4 alignment and dependency boundary

Added native `index` and `align-count` commands for explicit STAR/HISAT2 selection,
sequential argv subprocesses, default one worker, Slurm CPU-budget validation,
samtools auxiliary-thread accounting and explicit sort memory. Genome-only index
metadata binds FASTA/GTF hashes and STAR sizing parameters. Successful output is
published with Linux no-replace rename; errors and interruptions retain stage
logs and propagate failure without a completion marker. See
[P4_BACKENDS.md](P4_BACKENDS.md) for the full CLI and limits.

- `make -j2 check`: 36 cases / 222 assertions plus native/count/P4 integrations.
- Real STAR/HISAT2 × SE/PE × threads 1/2: all eight pass the independently stated
  P0 counts, assignment/ambiguity summaries, BAM checks and splice-junction CIGAR.
  The gate also passes under `unshare -Urn`, with networking disabled.
- ASan/UBSan full native gate passes; focused sanitizer gate rerun after the
  no-replace publication fix. Review fixes for private child PATH (staged wrapper
  interpreters) and exact version tokens pass their RED/GREEN regressions,
  scoped sanitizer checks and independent re-review. All eight real combinations
  passed again under disabled networking with ambient PATH=/usr/bin:/bin.
- `tools/toolchain.py stage` reconstructs the pinned alignment environment into
  a new prefix; all 48 installed package records match the lock.
- `tools/toolchain.py preflight` verifies full package locks, recipe/license
  hashes, executable versions, pinned R packages, DESeq2/apeglm and PNG output.
  It passes with network disabled using the newly staged alignment prefix.
  Five Python safety checks pass. R user libraries and auxiliary thread counts
  are constrained by preflight; compute never invokes installation.
- `vendor/toolchain/manifest.json` retains core package source URL/hash, binary
  hash, license, dependency and upstream recipe/build evidence. Binary packages
  are reproduced from locks; a source rebuild is not claimed. In particular the
  historical HISAT2 recipe's extra simde clone is unpinned. R post-link downloads
  mean environment staging remains online, distinct from verified offline compute.

The user confirms real data from the repository report for P6. The report names
GSE80336, a 36-sample cohort and a separate 35-sample analysis excluding C_28.
[P6_REPORT_DATA.md](P6_REPORT_DATA.md) records provenance, HTSeq-count historical
differences, and contrast-direction reconciliation needed before comparison.

### P5 workflow complete

`workflow-plan`, `workflow-local`, `workflow-submit`, and `workflow-task` now join
P1–P4 with immutable snapshots, hash-verified resume, shared index locks, failure
quarantine, cancellation and explicit Slurm dependencies/resources. Final native
gate passed (36 cases/222 assertions plus integration gates); real P0 alignment
and merge passed, with the genuine tiny-data R failure retained. Review fixes
for early Slurm job timing, separate parsable stdout, isolated pinned R preflight
and checksum-child cancellation passed focused regressions and final review.
See [P5_WORKFLOW.md](P5_WORKFLOW.md). Live Slurm execution remains untested.

### P6 reduced local check complete

The user explicitly requests this local 8 GB RAM / 8 CPU / 500 GB machine and
a 30-minute reduced check. Staged six samples × 50,000 prefix reads, with the
exact NC_000022.11 chromosome from checksum-verified report references (981 genes).
The local comparison uses two threads per tool, one R worker and two concurrent
local STAR alignment tasks after resuming a slower pilot. A 10k-read thinning
attempt failed normalization and was retained. One functional comparison completed
within the deadline; three repetitions are deferred. Retain HISAT2 as
default. Restricted reference/prefix sampling
cannot establish full-genome scientific or performance equivalence.

Published-count checks completed: 12 significant genes for full36 and 15 excluding
C_28, matching report totals, with optimizer warnings retained. The independent
six-sample DESeq2/apeglm numerical oracle passed after correcting incidental
vector-name attributes in its comparison helper; no production R change was needed.
Three source/staging checks also passed.

Both 50k-read workflows completed through real R and plots at 13:38:50 UTC,
27m15s within the 30-minute budget. Counts: HISAT2 877, STAR 985; Pearson 0.96209;
both zero DEGs at this depth. Exact sample/gene identity and merged/assigned count
totals agree. One interrupted/resumed comparison does not establish repeatability
or a clean speed ranking. Full results and limitations: [P6_REPORT_DATA.md](P6_REPORT_DATA.md).

### Next step

P7 migration/cleanup remains next; full-data repeated benchmarking is deferred.
Keep HISAT2
as default until the planned P6 comparison. Installation/build/test details and
remaining limits are in [LOCAL_DEVELOPMENT.md](LOCAL_DEVELOPMENT.md).

## Initial session record (superseded status)

## Serena and environment

Serena 1.5.3 responds with RNASeqDGE active, LSP backend, Bash language,
editing/interactive modes. `initial_instructions`, `get_current_config`, and
`get_symbols_overview(relative_path="run_pipeline.sh")` succeeded. Symbol lookup
returned the launcher variables; the previous Node PATH error did not recur.
This checks navigation, not every Serena operation.

Read the handoff, execution plan, migration rationale, legacy merger/metadata,
pipeline, R analysis, launcher, example config, and existing cluster smoke tests.
No repository AGENTS.md was found; the user-provided Context7 instructions apply.
This work uses source inspection and new standard-library scripts, with no new
third-party API implementation or dependency version pinning.

Local inventory (`tests/output/p0/preflight.json`, ignored generated evidence):

| Component | Observed |
| --- | --- |
| OS | x86_64 Linux, WSL2 kernel 6.18.33.2, glibc 2.43 |
| GCC / Make | 15.2.0 / 4.4.1; no native project build yet |
| Python | 3.13.13, `/home/dzk55/miniforge3/bin/python3` |
| R | 4.5.2 |
| Node | 24.21.0 |
| STAR, HISAT2, featureCounts, samtools | not found on current PATH |
| Snakemake, fasterq-dump, sbatch | not found on current PATH |
| pandas, biomart, loguru, Bio, snakemake Python modules | unavailable in selected interpreter |
| Required R packages | all 12 inventoried packages unavailable; see JSON |

These are observations about the current shell, not proof that no other machine
or environment has the tools. Cluster partition/QoS, scratch, memory limits, and
production sample mapping remain unknown. Version inventory is not a pinned
toolchain or an R execution preflight.

## Added artifacts and verified results

- `tools/p0_preflight.py`: read-only inventory with command timeouts and captured
  versions/failures. No installation, reference downloads, or scheduler calls.
- `tools/p0_legacy_probe.py`: records actual merger startup and isolated original
  awk/R expressions. Does not execute the workflow or query metadata services.
- `tests/fixtures/p0/`: four samples, five explicitly mapped technical runs,
  four genes, optional annotation, explicit comparison/settings, five shuffled
  count inputs, and a manually derived expected integer matrix. Fixture README
  explains each expected sum and keeps it distinct from legacy output.
- `tools/generate_p0_sequences.py`: deterministic 3 kb synthetic reference,
  exon/GTF records, 75-base SE reads and PE fragment, origins, reference hashes,
  and production-schema run/reference manifests. Includes plus/minus, splice
  junction, overlap, and zero-coverage cases. These reads are independent of
  the hand-authored count oracle and have not been aligned or counted.
- `tests/integration/check_p0_fixture.py`: **3 tests pass**. Checks independent
  arithmetic, complete gene sets, sample ordering, run/sample identities,
  technical-run aggregation, comparison/annotation settings, reference hashes,
  FASTQ lengths/qualities, and read sequences against recorded genomic origins.

Actual legacy evidence (`tests/output/p0/legacy.json`, separate from expected data):

1. `merge_transcripts.py --help` exits 1 at `import pandas` with
   `ModuleNotFoundError`. No merger matrix or scientific baseline was generated.
2. The original count-rule Python string's awk program emits a literal tab and
   exactly matches `gene_A.1<TAB>7<LF>` on a tiny input. This isolated check does
   not establish that Snakemake executes the entire rule successfully.
3. R evaluates the original condition expression: `control-2`, `bipolar-4`
   become `control`, `bipolar`. Alphabetical reference is **bipolar**, not control.
   The suspected regex failure was not reproduced on this R version.

Source inspection also shows the annotation inner join, missing mapping skip,
GEO sample intersection, and differing 0.05 result / 0.1 volcano thresholds.
Runtime gene loss, full DAG target behavior, fasterq-dump options, and statistical
contrast behavior have **not** yet been reproduced with pinned dependencies.

## Reproduction

Run from the repository root (Python stdlib; R/awk used by the legacy probe):

```sh
mkdir -p tests/output/p0
python3 tools/p0_preflight.py > tests/output/p0/preflight.json
python3 tools/p0_legacy_probe.py > tests/output/p0/legacy.json
python3 tests/integration/check_p0_fixture.py
git diff --check
```

All commands above were run. Inventory/probe commands succeeded in recording
evidence; this does not turn their recorded missing-tool/startup failures into
successful pipeline runs. Fixture tests and whitespace checks passed.

To deliberately regenerate only synthetic sequence inputs and their manifests:
`python3 tools/generate_p0_sequences.py` (also run during initial creation).
The hand-authored expected matrix is never regenerated by this command.

## Remaining gate and next work

P0 requires an executable original-versus-reviewed baseline, not just a fixture.
Before P1:

1. Stage a pinned legacy Python/Snakemake environment and reproduce merger
   gene/sample loss using cached synthetic annotation/mapping without network.
   Preserve original failure/output separately from any corrected baseline.
2. Validate synthetic FASTA/GTF/reads with staged HISAT2 2.2.1 and featureCounts;
   review explicit counting policy and observed assignment against read origins.
   Resolve tiny-index settings with current upstream docs/Context7 before use.
3. Finalize analysis key/value validation and covariate type/level encoding
   (the fixture currently covers categorical `condition` only).
4. Establish a pinned R/Bioconductor environment for the later P3 gate. The tiny
   arithmetic matrix is not a DESeq2/VST parity fixture; a suitable statistical
   fixture and corrected statistical baseline remain to be built.

P0 is intentionally reported as incomplete. Missing dependencies prevent a
scientific baseline claim; no legacy output has been promoted to golden data.
