# P5 local and Slurm workflow

Execute after P4 review fixes. User explicitly authorizes continuing P5 and P6,
on the current dirty main without commits. Use P1 bundle validation, P2 merge,
P3 offline R and P4 native index/align-count. HISAT2 default is retained.

## Task 1: reproducible native workflow and submission boundary

Add native workflow commands and a small local Make entry point. Own
src/workflow.cpp, include/rnaseq/workflow.hpp, src/main.cpp, Makefile,
tests/integration/check_p5_workflow.sh, optional tests/integration/check_p5_real.py,
docs/P5_WORKFLOW.md and new example config/manifest files. Existing P1–P4 work
must be preserved. No subagents or commits.

Required behavior:

- Inputs are explicit canonical samples/runs/references/analysis/contrasts,
  explicit gene universe and optional annotation, tool bin dir, pinned Rscript,
  backend default hisat2, threads/workers default1, STAR sizing, run directory,
  optional shared index-cache directory. No network, installs, automatic study
  metadata or implicit FASTQ/reference download. Validate all inputs and resource
  budgets before work. R workers and native threads obey Slurm allocation and
  auxiliary BLAS/OpenMP caps; stages are serial locally regardless of make -j.
- `plan` or equivalent prints DAG, resolved paths and requested CPU/memory knobs
  without creating outputs, launching analysis tools, downloading or submitting.
- Prepare an immutable per-run input snapshot with resolved data paths, SHA-256
  for all inputs, executable/script identity and relevant package locks; keep
  metadata copies. Large FASTQ/reference contents remain external and are hashed
  and rechecked rather than copied. No overwrite of global config. On attempted
  resume with changed source/config/tool hashes fail clearly and require a new
  run directory; this is safe invalidation, never reuse stale downstream output.
- Stages: shared/reference index -> each technical run alignment/counting ->
  canonical sample merge -> offline R. Build P2 inputs.tsv from featureCounts
  exact recorded BAM column names, including retained stage path in that header.
  Preserve all raw logs/summaries. Every successful stage has a completion
  manifest of stage identity and output hashes, written last via atomic publish.
  Resume validates these hashes, not file existence or mtime. Partial/corrupt
  output must not be silently reused; quarantine invalid stage directories while
  preserving diagnostics, rerun failed stage and dependent stages.
- Lock the run to prevent concurrent local drivers/preparations; support Slurm
  per-run task concurrency with individual stage locks and consistent snapshots.
  Lock reference index cache by complete key (reference, backend/tool version,
  all index params); publish only complete cache entries. If cache support would
  compromise correctness, report a precise design conflict before omitting it.
- Cancellation must terminate active descendants and leave no success manifest.
  Existing P4 child process-group behavior must remain correct when the workflow
  calls native subprocesses (avoid leaving grandchildren in separate groups).
- Slurm: explicit submit command uses argv `sbatch --parsable`, creates index job
  -> alignment/counting array with bounded `%concurrency` -> merge/R afterok job.
  Job scripts select array task from immutable run mapping. Explicit CPU/memory/
  time options, optional partition; no nested submission. Record accepted job IDs
  immediately and refuse accidental duplicate submit. On partial submission
  failure cancel already accepted jobs using argv scancel and retain diagnostics.
  A generated plan is not reported as a real cluster execution.
- Local Make target passes a single config/path to the native workflow; shell
  quoting must preserve literal paths, no user strings interpolated into shell
  expressions. Keep basic `make check` offline and free of R/Python/alignment
  requirements.

Test with public CLI and lightweight tool/sbatch/scancel/R stubs: successful DAG,
no-op verified resume, changed source and config/tool invalidation, corrupt and
partial output recovery, failure/restart, concurrent run/cache protection,
cancellation including descendants, read-only plan, Slurm exact dependency/array
resources, no duplicate submit and partial submission rollback. Add a small real
local integration where meaningful (P0 read fixture is too small for DESeq2;
do not invent R success for it). Reuse P3's independently verified real R suite
for statistical behavior and state end-to-end limitations explicitly.

Use Context7 Slurm docs: --dependency=afterok:ARRAY_JOB_ID waits for the entire
array to succeed, --array=0-N%K bounds concurrent tasks, --parsable prints accepted
job identity. No scheduler is currently installed locally, so test submission
protocol with stubs and document live-cluster validation as pending.

Write report `.superpowers/sdd/P5_IMPLEMENTATION_PLAN/task-1-report.md` with exact
commands/outcomes, scoped review concerns and limitations. Controller owns
progress/handoff/Serena and P6 data staging; avoid touching those files.
