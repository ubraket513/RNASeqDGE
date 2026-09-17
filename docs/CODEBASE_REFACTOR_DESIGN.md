# Maintainability refactor design

Status: proposed for review, 2026-09-17. This is a behavior-preserving refactor
of the existing C++20/R/Bash workflow in the current main checkout.

## Findings

- `Makefile` repeats the same compilation recipe for individual source and test
  files. Its C++20 requirement disagrees with the README's C++17 requirement.
- `src/workflow.cpp` combines configuration parsing, resource allocation,
  immutable snapshots, stage publication, profiling, local execution and Slurm
  submission. Its 356 physical lines hide substantial complexity: 22 exceed
  180 characters, with many statements and branches packed onto single lines.
- Preparation/checking R scripts contain 70 lines longer than 180 characters.
  Nested source-path expressions and dense control flow make review difficult.
- Current usage guides and historical implementation plans share one docs
  directory. Names such as P4/P5 explain development chronology, not user tasks.
- There is no checked-in editor/formatter policy or automated CI entry point.
- Existing independent scientific oracles and process/cancellation tests are
  valuable constraints. They should remain independent of production helpers.
- Workflow identity hashes the R entry script but does not discover sourced R
  modules. Splitting production R blindly would weaken resume invalidation.

## Approaches

1. Formatting and documentation only: smallest risk, but leaves workflow
   responsibilities coupled.
2. Incremental internal refactor with stable public paths: recommended. Improves
   module boundaries and build consistency without forcing new user interfaces.
3. Wholesale build-system/package/CLI replacement: more migration effort and
   new failure modes without evidence that the current project needs it.

## Recommended design

Keep `src/`, `include/rnaseq/`, `tests/unit/`, `tests/integration/`,
`tests/fixtures/`, `tools/`, `config/`, `examples/`, `third_party/` and `vendor/`.
Keep the existing launcher, native commands and R CLI paths. These are already
reasonable conventions for a small native scientific application.

Introduce private workflow implementation files under `src/workflow/`:

- `config`: configuration validation, path resolution and resource decisions.
- `state`: locking, snapshots, inventories, verified stage publication and
  diagnostic profiling.
- `execution`: index/alignment/merge/R stage coordination and Slurm submission.

Keep `src/workflow.cpp` as command dispatch and the public header as the narrow
entry interface. Use concrete internal types and existing process supervision;
do not add a generic plugin framework or replace fork-based concurrency.
Extract in small steps, preserving validation order and source identities.

Simplify Make with explicit source lists and shared compile rules, separate
native and test object paths, and dependency files. Preserve compiler/linker
overrides, offline `check`, sanitizer builds and alternate `BUILD` directories.
Keep GNU Make; changing build systems is not needed for this cleanup.

Add a small editor/format policy. Expand compressed C++ and R statements into
readable control flow and name nontrivial intermediate values. Consolidate
duplicated helper code only where behavior and ownership are genuinely shared.
Keep scientific formulas, thresholds, filters and independent test calculations
unchanged. Keep the production R script self-contained in this pass so its
complete code remains covered by the existing snapshot hash.

Give `docs/` a short index separating operating/development guides, validation
evidence and historical plans. Move historical planning material into
`docs/history/`, update references, and preserve scientific reports, sample
results, fixtures, package provenance and license notices. Do not delete data
merely because it is not part of the executable.

Add concise contributor instructions and a minimal native CI gate using
`make check`; installed-tool/R tests stay opt-in with documented commands.
CI must not fetch research data or install Python. Add no production dependency.

## Compatibility and verification

Preserve commands, TSV schemas, output names, numerical tolerances, atomic
publication, failure retention, verified resume, resource semantics, signal
propagation and scheduler dependencies. Rebuilt executable/source hashes will
invalidate old workflow snapshots as they already do; never bypass this check.
Leave the existing README/license edits intact and do not commit automatically.

Run the native baseline before edits, then targeted tests after each extraction
and the complete native gate at completion. Exercise alternate build directories
and sanitizer checks for touched process/state code. For R formatting, compare
parsed expressions before/after when applicable, run helper checks and preserve
the independent DESeq2 parity gate. Rerun a small installed-tool workflow when
execution boundaries change. Check documentation links and whitespace.

Success means clearer ownership and less duplication with the same validated
behavior. It does not mean fewer lines at the expense of readability, higher
performance without measurements, or deleting reproducibility evidence.
