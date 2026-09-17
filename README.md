# RNASeqDGE

A native C++ RNA-seq workflow with offline R/DESeq2 analysis. The compute runtime
uses no Python interpreter. It aligns staged FASTQ reads with HISAT2 or STAR,
counts with featureCounts, merges canonical samples, and runs DESeq2/apeglm.
The research report studies bipolar disorder using GSE80336 / PRJNA318642.

## Build and prepare

Linux, a C++17 compiler, Make, Bash and standard Unix utilities are required.
Build and run the offline native checks:

```sh
make -j2
make check
```

Prepare pinned packages once, using a standalone `mamba` or `micromamba` executable
(pass `--mamba /path/to/micromamba` if needed). Commands require absent destination
prefixes; do not rerun them over an existing environment.

```sh
mkdir -p .deps tests/output
bash tools/toolchain.sh stage runtime-r --prefix .deps/runtime-r
export PATH="$PWD/.deps/runtime-r/bin:$PATH"
bash tools/toolchain.sh stage alignment --prefix .deps/alignment
bash tools/stage_native_runtime.sh .deps/alignment .deps/runtime-tools
Rscript --vanilla tools/toolchain.R preflight \
  --alignment .deps/runtime-tools --r-prefix .deps/runtime-r \
  --out tests/output/runtime-preflight
```

The alignment **source staging prefix** contains upstream wrapper dependencies,
including Python. The extracted `.deps/runtime-tools` payload contains native
executables and their private libraries, plus STAR's Bash dispatcher; compute
jobs use this payload and `.deps/runtime-r`, not the source staging prefix.
The R lock contains no Python package. System Python is neither used nor removed.
Historical P0 locks remain provenance and are not the compute environment.

HISAT2 invokes its native small-index binaries directly with the upstream wrapper
marker. This interface accepts uncompressed FASTQ and small `.ht2` indexes;
compressed reads and large HISAT2 indexes are not supported by this workflow.

## Run

Create validated sample/run/reference/analysis/contrast manifests and a two-column
`key`/`value` workflow TSV. See [workflow configuration](docs/P5_WORKFLOW.md)
and [the example](examples/workflow-p0.tsv). Paths resolve relative to their
configuration or manifest. Data and references must already be staged.

```sh
bash run_pipeline.sh plan study-workflow.tsv
bash run_pipeline.sh local study-workflow.tsv
# Explicitly submit on a configured Slurm host:
bash run_pipeline.sh submit study-workflow.tsv
```

Planning validates inputs without running analysis. Runs retain immutable input
identities, logs, counts, statistical tables, plots, and per-stage profiles.
Reusing the same configuration resumes only verified stages. Changed inputs or
software require a new run directory. Failed outputs remain available for diagnosis.
The tiny P0 example intentionally reaches a too-few-genes R failure; use a real
study for successful differential expression.

Local sample concurrency is bounded by `local_cpus`, `threads`, `local_mem_mb`
and `local_job_mem_mb`, with optional `local_jobs`. Without both memory estimates,
alignment is serial. Estimates are reservations, not enforced memory limits.
R `workers` is separate; BLAS/OpenMP helper threads are capped at one. Compatible
contrasts reuse a fitted DESeq2 model within an invocation. DESeq2 remains in R
with its compiled numerical routines, preserving the validated statistical method.

## Evidence and limits

[Progress](docs/IMPLEMENTATION_PROGRESS.md), [local development](docs/LOCAL_DEVELOPMENT.md),
[Python-free runtime validation](docs/P7_RUNTIME.md), and
[report-data validation](docs/P6_REPORT_DATA.md) describe the checks and evidence.
The six-sample, chromosome-22 comparison is a reduced functional check, not a
whole-genome performance benchmark or reproduction of the report's biological conclusions.
Live Slurm execution still needs validation on a cluster.

The former Python/Snakemake entry points are retired; see
[legacy recovery](docs/LEGACY_RECOVERY.md). Report sources and historical evidence
are retained.
