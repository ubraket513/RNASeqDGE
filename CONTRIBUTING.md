# Contributing

Use GCC 12+ (C++20), GNU Make, Bash and standard Linux utilities. Native tests
run offline and do not require R, Python or installed alignment tools.

```sh
make -j2 check
# Independent build directory; useful for compiler/flag changes:
make -j2 BUILD=build/debug CXXFLAGS='-O0 -g' check
# Address/undefined-behavior instrumentation:
make -j2 sanitize
```

## Code layout

| Path | Responsibility |
| --- | --- |
| `src/`, `include/rnaseq/` | Native implementations and public interfaces |
| `src/workflow/` | Private configuration, stage state and execution modules |
| `run_deg_analysis_offline.R` | Self-contained statistical entry point |
| `tools/` | Preparation, runtime packaging and diagnostic helpers |
| `tests/unit/`, `tests/integration/` | Native unit tests and boundary/scientific gates |
| `tests/fixtures/` | Small deterministic fixtures and independent expectations |
| `config/`, `examples/` | Pinned dependency provenance and workflow examples |
| `docs/`, `docs/history/` | Current guides/evidence and historical plans |
| `report/`, `sample/` | Research report and retained example results |
| `third_party/`, `vendor/` | Upstream code, recipe provenance and license notices |

Keep public commands, script locations, TSV contracts and scientific outputs
stable. Prefer small concrete modules over generic frameworks. Internal workflow
headers stay under `src/workflow/`; avoid exposing implementation details through
`include/rnaseq/`.

## Style and change discipline

Follow `.editorconfig`: four-space C++ indentation, two-space R/shell indentation,
UTF-8/LF, and tabs in Make recipes. `.clang-format` defines C++ formatting. Format
only owned source files, never vendored code. Expand branches and loops so reviewers
can follow control flow; line count is not a measure of simplicity.

Keep formatting, behavior changes and scientific changes reviewable separately.
For mechanical R formatting, compare `parse(..., keep.source = FALSE)` before and
after. Do not use shared production calculations as the independent test oracle.

A workflow snapshots executable and R-script hashes. Rebuilding or editing source
while a real workflow is active deliberately invalidates it. Use a new run directory
after software changes. The R entry point stays self-contained so its hash covers
its implementation. Local memory settings are reservations, not enforced RSS caps.

## Verification

Run `make check` and `git diff --check` for native changes. Add focused regressions
for new behavior; use existing coverage for mechanical extraction. Run sanitizer
checks when changing process/state code. The native CI executes the same offline
gate on Ubuntu; installed-tool and scientific tests are separate opt-in gates.

For R changes, use the pinned runtime and relevant independent oracle:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p7_fit_cache.R
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p6_staging.R
```

Real-tool checks and the complete statistical suite are documented in
[local development](docs/LOCAL_DEVELOPMENT.md). Schedule memory-heavy tests
sequentially on an 8 GB host. Preserve scientific tolerances and retain logs for
failed runs. Do not infer whole-study or live-Slurm validation from stub tests.

Project code is MIT licensed; preserve [third-party notices](README.md#license).
Do not commit dependency prefixes, generated run directories or raw research data.
