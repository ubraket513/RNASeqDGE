# P5 native workflow

`rnaseq workflow-local CONFIG.tsv` executes the DAG: index, bounded parallel technical
run alignment/counting, canonical sample merge, then the offline R command.
It never fetches data, installs tools, changes global configuration, or discovers
study metadata. P1–P4 scientific interfaces remain unchanged. HISAT2 is the default;
STAR is selected explicitly. `make -j` does not control workflow concurrency.

```sh
make
build/rnaseq workflow-plan examples/workflow-p0.tsv
WORKFLOW_CONFIG=examples/workflow-p0.tsv make -j4 workflow
```

The P0 example has real reads but too few genes for DESeq2: it demonstrates genuine
alignment/merge followed by a retained R failure, **not** an end-to-end successful
scientific analysis. Supply an adequately sized study for a successful R stage.

The strict configuration schema has exactly `key` and `value` TSV columns.
Required paths are `samples`, `runs`, `references`, `analysis`, `contrasts`,
`genes`, `bin_dir`, `rscript`, `r_script`, `tool_lock`, `r_lock`, and `run_dir`.
`annotation` and `index_cache` are optional paths. Paths resolve relative to their
configuration file; embedded run/reference paths resolve relative to their own
manifests. All scientific manifests use the existing P1/P2/P3 contracts. Large
FASTQ and reference files must already exist locally and FASTQs must be uncompressed.
The tool locks are explicit package locks or the native bundle's `SHA256SUMS`
inventory retained as provenance. Executable
versions are checked by P4; P5 explicitly checks R 4.5.3 and DESeq2 1.50.2, apeglm 1.32.0,
BiocParallel 1.44.0 and pheatmap 1.0.13 before its analysis invocation.
R_HOME, R_LIBS, LD_LIBRARY_PATH and LD_PRELOAD overrides are cleared for R;
user/site R libraries point to /dev/null, keeping packages in the selected prefix. Locks alone do
not audit every installed dynamic library or every installed R package file.

Optional settings are `backend` (`hisat2`/`star`), `threads` (1), `workers` (1),
`star_sa_bases` (14, range 1–14), and `star_chr_bits` (18, range 1–18). Native threads
and R workers must each fit `SLURM_CPUS_PER_TASK` and the configured Slurm CPUs.
BLAS/OpenMP auxiliary worker environment variables are capped at one. Alignment
and R execute in separate phases, so their worker requests are not summed. P4 bounds
samtools sort memory per thread; this is not a total RSS guarantee. Scheduler
memory and time settings are allocation requests rather than local process limits.

Local scheduling accepts `local_jobs` (maximum simultaneous alignments),
`local_cpus` (total CPU budget), `local_mem_mb` (total reservation budget), and
`local_job_mem_mb` (estimated peak memory of one complete alignment/counting job).
CPU defaults to process affinity and is capped by `SLURM_CPUS_PER_TASK`; an explicit
`local_cpus` can lower it. R `workers` and alignment `threads` must each fit that
resolved CPU budget. Slurm `SLURM_MEM_PER_NODE` or `SLURM_MEM_PER_CPU` also caps
available reservations when present. Positive integers are required.

With both total and per-job memory estimates, resolved concurrency is the minimum
of the requested jobs (default CPU-derived), CPU budget divided by `threads`, and
memory budget divided by the per-job estimate. Thus concurrent alignment threads
and reservations cannot exceed their configured budgets. An oversized single-job
reservation is rejected. Missing either memory estimate gives a clearly reported
serial fallback. These estimates are user-supplied reservations, **not an RSS
limit**: profile representative data before enabling parallelism, and account for
index size, sorting, counting, and filesystem cache needs. Index and R stages remain
single DAG phases; the alignment reservation does not estimate R worker memory.

For example, `threads=4`, `local_cpus=16`, `local_jobs=4`, `local_mem_mb=32000`,
`local_job_mem_mb=12000` resolves to two simultaneous alignments, reserving 8 CPUs
and 24000 MB. A smaller Slurm CPU allocation can reduce that further.

Submission resource validation uses the requested `slurm_cpus` and
`slurm_mem_mb`, independently of the login host's affinity or local budgets.
`workflow-plan` reports `resource_scope=slurm_requested` when `slurm_cpus` is
configured; otherwise it reports the local allocation. `workflow-local` and
executing `workflow-task` commands still enforce actual affinity and allocation
limits before work starts. Planning a Slurm submission therefore does not prove
that the same configuration can run locally.

`workflow-plan` validates and hashes inputs, prints resolved paths, requested
resource settings and the DAG, and creates no output. Hashing invokes the existing
native SHA-256 helper (`sha256sum`); planning runs no analysis/version command,
network operation, or scheduler command.

## Snapshot, completion, and resume

The run is protected by a nonblocking advisory `flock`. A local driver or submitter
requires an exclusive lock; scheduler tasks wait for a shared run lock while preparation/submission completes,
then take individual
stage locks. Contention fails with `busy lock` and is safe to retry. This requires
a filesystem with coherent cross-node `flock` behavior for live cluster use.

`snapshot/` contains copies of metadata and package locks, resolved configuration,
resolved reference/run tables, and SHA-256 identities for every declared source,
the native executable, R script/Rscript, and regular executable-directory files
(including dereferenced wrapper companions). Large reads/references stay external.
Snapshot publication is an atomic directory rename; its hash manifest is written
last. Subsequent commands rehash the sources and snapshot. Changed source,
configuration, tool or snapshot content requires a new run directory. Sources are
also rehashed before publishing each successful stage. Users must keep input
files stable during a run; hashing is not a transactional filesystem snapshot.

The checksum subprocess uses Linux parent-death signalling so cancellation during
large-file hashing stops the checksum reader as well.

Each successful output has `STAGE.tsv`, containing its dependency identity and
hashes of every output file, including raw logs, summaries and provenance.
`generation.txt` changes whenever a stage is rebuilt, forcing downstream rebuilds
even if regenerated scientific values happen to be identical. Verified stages
are reused without executing tools. Missing/changed/extra files invalidate the
stage; it is renamed to `.invalid-PID[-N]` and rebuilt. P4 and R also preserve their
own failed staging directories. Failed R invocation logs are retained at run level.
No workflow success manifest is published after a failed or cancelled stage.
P4 `COMPLETE` is a lower-level publication marker; P5 trusts only `STAGE.tsv`.

P2 `merge/inputs.tsv` reads the exact single BAM column from each featureCounts
header, including its historical P4 staging path. It does not invent a final-path
BAM column. Both technical counts and raw featureCounts summaries remain available.

An optional shared index cache is keyed by reference hashes, genome-only policy,
backend, executable/lock identities, native executable identity, thread count,
and all STAR index parameters. Keys deliberately over-invalidate across native
rebuilds or unrelated selected-tool changes. Exclusive per-key locks protect
build/publication; alignments hold shared locks while using an index. Partial or
corrupt entries are quarantined. Cache checks verify file hashes, not mtimes.
A competing cache builder receives `busy lock`; no unsafe fallback index is used.

Local alignment jobs run in forked workers, each with the existing P4 process
supervisor. No C++ threads compete for process-global signal handlers. Stage/cache
locks and verified resume stay in effect. A failed worker stops its siblings and
prevents merge/R. The parent traps SIGINT/SIGTERM, stops and kills the descendant
tree (including tools that changed process groups), and reaps it using Linux
subreaper support before returning 130/143. Cancellation during index/R retains the
existing synchronous supervisor. These Linux facilities match the existing
`/proc/self/exe` runtime requirement.

Each stage attempt adds a TSV under `run_dir/profiles/`, outside hashed scientific
outputs. It distinguishes `executed`, `reused`, and ordinary `failed` attempts,
separates stage validation, execution, source revalidation, and publication hashing,
and records wall time, CPU deltas, requested reservations and resolved concurrency.
RSS columns are **process-lifetime high-water marks** from Linux `getrusage` for
self and waited-for children; they are not a summed concurrent peak or stage RSS
delta. CPU values include waited-for checksum/tool descendants. Snapshot/preflight
hashing and index validation outside the stage function are not attributed to these
stage rows. Abruptly killed workers may not publish a profile row; retained logs
and absent stage markers remain the failure evidence. Profiles are diagnostic and
do not affect `STAGE.tsv` integrity or scientific tables.

## Explicit Slurm submission

Add `slurm_bin_dir`, `slurm_cpus`, `slurm_mem_mb` (positive integer MB), and
`slurm_time` to the config, plus optional `slurm_partition` and
`slurm_concurrency` (default 1). The binary directory must contain executable
`sbatch` and `scancel`. Then explicitly run:

```sh
build/rnaseq workflow-submit study-workflow.tsv
```

This prepares the snapshot and submits three argv-only commands using
`sbatch --parsable`: index → alignment/counting array `0-(N-1)%K` → merge/R.
The array depends on `afterok:INDEX_ID`; merge/R depends on `afterok:ARRAY_ID`,
which requires every array task to succeed. CPU, memory, time and optional partition
are passed explicitly to every job. Scripts safely quote literal paths and call
`workflow-task CONFIG index|array|finish`; array tasks select rows from the validated,
immutable snapshot mapping using `SLURM_ARRAY_TASK_ID`. No nested submission is
allowed from a process with `SLURM_JOB_ID`.

`submission/` retains scripts, exact argv, separate stdout/stderr submission logs and immediately recorded
accepted IDs. Any attempted submission prevents an accidental second submission
or a local driver on that run. A partial submission failure calls `scancel` for
all known accepted IDs and retains cancellation outcomes. Inspect scheduler state
and diagnostics before preparing a new run after any submission interruption;
a scheduler accepting a job before an interrupted client receives its ID is an
inherently ambiguous submission boundary. Generated scripts and accepted IDs are
not evidence that a cluster job actually executed.

The Slurm protocol follows [job-array dependencies](https://slurm.schedmd.com/job_array.html)
and [sbatch options](https://slurm.schedmd.com/sbatch.html), checked with Context7.
Live multi-node filesystem locks, scheduler execution, cluster routing and memory
sizing remain site-specific validation tasks. No live scheduler is installed in
this workspace; only the submission protocol and generated task scripts have been
executed with stubs/local task environments.

## Verification

`make check` remains offline and needs no R/Python/alignment installation. Its P5
shell stubs test CPU/memory resolution, actual bounded overlap, failed-sibling and
signal cancellation with separately grouped descendants, deterministic merge,
profiling, the DAG, read-only planning, resource rejection, no-op hash-verified
resume, source/config/tool invalidation, corruption and partial-output recovery,
failed R restart, cache reuse/contention, run locking, descendant cancellation,
Slurm resources/dependencies/array mapping, duplicate-submission rejection and
partial-submission rollback. Stub R output tests orchestration only.

```sh
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p5_real.R
```

The separate installed-tool gate uses real HISAT2/samtools/featureCounts on five
P0 technical runs and checks independently expected counts after their merge into
four canonical samples. It requires the actual too-small-data R failure and absence
of an R-stage success manifest. P3's independently verified full statistical suite
supplies the DESeq2/apeglm correctness evidence. A successful full raw-read-to-DE
study and a live Slurm run remain unvalidated by P5's tiny fixture. Subsequently,
the P6 reduced six-sample chr22 check completed through real R for both backends;
see [P6_REPORT_DATA.md](P6_REPORT_DATA.md). This does not validate a full-genome study.
