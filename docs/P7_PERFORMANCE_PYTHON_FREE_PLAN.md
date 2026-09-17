# Performance and Python-free runtime implementation plan

User approved improving scheduling/profiling/fit reuse while retaining DESeq2,
and explicitly requires eliminating the Python interpreter from compute runtime.
Work on the existing main checkout; no commits or unrelated cleanup. Current
baseline HEAD is 2f28776; preserve existing documentation/Serena changes.

## Design and constraints

Keep native C++ orchestration/count handling and pinned R/DESeq2/apeglm scientific
implementation. Introduce bounded process concurrency locally; never call the
signal-owning process supervisor from competing C++ threads. Preserve shared
index/stage locks, immutable source identity, verified resume and descendant
cancellation. Keep BLAS/OpenMP caps to prevent nested oversubscription.

Remove all repository-owned Python and Snakemake execution paths, replacing
required preparation/checks with shell/R/native equivalents. Invoke HISAT2's
native aligner directly after equivalence checks. Provide a compute dependency
environment without a Python interpreter, distinct from environment creation.
Retain licensed vendor provenance and historical report text as historical data.
Do not silently remove capability or weaken independent scientific gates.

## Tasks

1. Native local scheduler/profiling: own src/workflow.cpp, process supervisor
   extensions if needed, associated native tests/Makefile. New local resource
   settings express total CPU, total memory, per-job memory and max parallel jobs.
   Bound concurrency by both budgets; absent explicit memory estimates retain a
   safe serial fallback, clearly reported by plan. Stage wall/CPU/RSS and hash
   validation timings must distinguish reused from executed stages. Regression
   checks exercise overlap, budget rejection, failure sibling cancellation,
   signals, resume, and deterministic merge. No alignment scientific changes.
2. R fit reuse/profiling: own run_deg_analysis_offline.R and R tests plus benchmark
   helper. Cache only compatible models (identical design matrices/levels/counts);
   retain independent fits for reparameterized coefficients. Add auditable timings
   outside scientific tables. Preserve raw/shrunken values, NA masks, coefficients,
   filtering, plots, and defaults. Test cached versus uncached independently and
   benchmark workers1/2/4 on the published full36 data, respecting observed available
   memory and a bounded runtime; record unsupported/resource-limited settings.
3. Repository helper migration: port required tools/*.py and tests/integration/*.py
   to R/shell/native with meaningful existing checks preserved. Remove obsolete
   legacy-specific probes with documented historical Git recovery. Do not remove
   metadata.py/merge_transcripts.py/pipeline.smk/launcher until root switches docs
   and entry points. No modifications to workflow.cpp or production R script.
4. Controller integration: Python-free HISAT2 adapter, runtime dependency staging,
   native launcher/default documentation, removal of obsolete Python/Snakemake
   entry points/dependency requirements, consolidated real-tool verification and
   migration recovery notes. Update gates/progress/Serena with actual evidence.

## Verification and review

Native offline make check and targeted concurrency/signal sanitizers as relevant;
independent R parity including compatible-model cache behavior; real alignment
oracle with direct HISAT2 binary; Python-free compute smoke with executable
availability restrictions; meaningful replacement helper tests; fresh native
build and opt-in real-data stage. Profile before recommending worker settings.
Independent task reviews and final integration review; retain diagnostics and
report limitations without treating a prepared configuration as executed evidence.
