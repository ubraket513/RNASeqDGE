# Implementation progress

Updated 2026-09-17. Baseline HEAD remains
`9f181e2310196fc22ce1159e530b320c48b4592c`.

## Status

**P0 local gate passed; P1 native foundation implemented, full manifest semantic
validation pending. P2–P7 not started.** HISAT2 remains the default. Dependencies
were installed after explicit user authorization. No production pipeline files
were changed, jobs submitted, commits created, or original-study statistical
parity claimed.
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

P1 is not claimed complete: `validate-tsv` is syntax-only for run/reference/
analysis/contrast manifests. Their semantic/cross-file validators remain to be
implemented, followed by P2 streaming count adapters and merge. The current
table loader materializes inputs in memory and is not a large-count performance
implementation. R integer limits/model validation and statistical parity remain
P3. CPU STAR is installed but has not yet been validated as an alignment backend.

### Next step

Finish P1 semantic validation against `DATA_CONTRACT_V1.md`: run/sample mapping,
layout/strand consistency, paths/reference roles and hashes, declared covariates,
and contrast levels. Then implement P2 count adapters/aggregation against both
the hand-calculated oracle and the reviewed P0 count baseline. Keep HISAT2 as
default until the planned P6 comparison. Installation/build/test details and
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
