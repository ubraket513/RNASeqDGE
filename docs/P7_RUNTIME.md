# Python-free runtime and performance changes

The supported compute environment consists of `.deps/runtime-r` (150 pinned
packages) and `.deps/runtime-tools` (native alignment/counting payload). No Python
interpreter or libpython is present in either prefix. System Python and historical
source environments are not uninstalled. Preparation can use a source package
prefix containing upstream wrappers; it is excluded from compute configuration.

Use the setup and launcher commands in the root README. Audit a staged runtime:

```sh
bash tools/check_python_free_runtime.sh .deps/runtime-r .deps/runtime-tools
.deps/runtime-r/bin/Rscript --vanilla tools/toolchain.R preflight \
  --alignment .deps/runtime-tools --r-prefix .deps/runtime-r \
  --out tests/output/runtime-preflight-new
```

The native bundle records copied files and SHA256 hashes, original package
provenance, and shared-library dependencies. It requires the recorded Linux host
ABI libraries. It is not a portable static executable distribution.

## Verification on 2026-09-17

- Runtime inventory/library checks and pinned-package/R preflight passed.
- Eight real-tool cases passed: HISAT2 and STAR, single/paired reads, one/two
  threads; exact counts, assignment summaries, BAM integrity and junction CIGAR.
  Evidence: `tests/output/p7-native-real-unsandboxed/verified.json`.
- Six 50,000-read samples from the report study, chromosome 22 reference, completed
  indexing, two simultaneous alignment/counting jobs (two threads each), merge,
  and real R analysis using the new launcher and Python-free prefixes.
  PATH contained only selected Unix utilities and those prefixes; neither legacy
  dependency prefix nor a Python executable was on PATH.
  Evidence: `tests/output/p7-python-free-workflow`, its sibling configuration/log,
  and `tests/output/p7-runtime-preflight`.
  Merged counts matched P6 byte-for-byte; all nine scientific TSVs passed strict
  numeric/identity/NA-mask parity against the previous HISAT2 analysis.
- featureCounts 2.0.6 aborted inside the execution sandbox, including the original
  source installation. The identical command succeeded outside it; the real-tool
  and workflow checks therefore ran outside that sandbox. No tool version changed.

## Performance and statistical evidence

Local alignments use bounded process concurrency with explicit memory estimates;
see P5_WORKFLOW.md. Profile TSVs distinguish execution from reuse. Reported native
RSS is a process high-water mark, not aggregate concurrent memory.

R caches fits only for identical validated model data and design matrices within
one invocation. Reverse baselines require separate fits. Independent oracle checks
and comparison to pre-change outputs passed for all 27 scientific tables in a
three-contrast fixture. Numerical tolerance remains `1e-10 + 1e-7 * abs(reference)`.

The full36 published-count benchmark with one R worker completed in 92 seconds,
with sampled aggregate peak RSS 1,291,372 KiB. All nine scientific tables matched
its reference. Two/four-worker runs were explicitly skipped because available
headroom failed the conservative memory guard. These observations do not establish
multicore R speedup. Re-run `tools/benchmark_r_workers.sh` when adequate memory is
available; its total budget is 30 minutes and it records resource-limited modes.

The raw-read run is a reduced functional check, not a whole-genome benchmark.
Live Slurm execution remains untested on this local machine.

## Helper migration

Repository-owned Python and Snakemake files are removed. R/shell replacements
cover preparation, provenance, fixture generation and integration checks. Retained
GEO cohort files, thinned reads and local-stage reads/reference/annotation matched
their original outputs byte-for-byte; comparison summaries matched numerically.
The local restaging process wrote its complete outputs but subsequently exited
nonzero after its source was edited while it was running. Its output parity is
evidence of the generated data, not a clean process exit. Final source parsing and
focused sampling-provenance regression checks passed. Other helper checks and the
independent complete compute workflow passed normally.
