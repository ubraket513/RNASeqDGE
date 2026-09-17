# C++ migration decisions

Execution plan: [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md).
Next-session entry point: [SESSION_HANDOFF.md](SESSION_HANDOFF.md).
The execution plan supersedes the preliminary implementation order below and
includes the subsequently accepted CPU STAR evaluation. Implementation has not started.

Status: direction accepted; dependency research completed; implementation pending.
Reviewed 2026-09-17. The first milestone replaces metadata/count-processing
utilities and workflow deployment while retaining a pinned R/DESeq2 backend.

## Reference implementation

Inspected [ubraket513/AntRepCLA](https://github.com/ubraket513/AntRepCLA/tree/806e2b0bd9d08f06c740c1e3bae80e177a3986ce)
at commit `806e2b0bd9d08f06c740c1e3bae80e177a3986ce`.
Relevant files are `src/igblast.cpp`, `src/types.hpp`, `src/export.cpp`,
`Makefile`, `tools/verify_output.sh`, and `docs/DECISIONS.md`.

Local verification with GCC 15.2.0: `make NPROC=2 check` passed all 36 tests
(108 assertions) and the frozen-output comparisons. Separate invocations with
explicit `--workers 1`, `2`, `4`, and `8` also matched the golden lineage tables
and summary statistics. This verifies the reference application, not RNASeqDGE.

Its useful model is a small C++ application with vendored headers, a direct
Make build, separate external sequencing executables, and frozen scientific
outputs. RNASeqDGE should adopt this packaging and validation approach without
assuming that the antibody-analysis performance measurements transfer to RNA-seq.

## Library choices

| Component | Observed in AntRepCLA | RNASeqDGE decision |
| --- | --- | --- |
| TSV/CSV parser | vincentlaucsb/csv-parser 5.3.0, amalgamated `csv.hpp` | Adopt a pinned upstream header for sample sheets and count tables. |
| Hash maps/sets | martinus/unordered_dense 5.0.1 | Candidate for gene-ID indexing; use only if profiling justifies replacing standard containers. |
| Unit tests | doctest header reporting 2.5.0 | Adopt a pinned header for parser, join, and validation tests; test binary only. |
| Parallel computation | OpenMP via GCC | Add only to measured parallel work; explicitly budget parser threads too. |
| CLI/logging | Project-owned code | Keep small initially; avoid adding libraries without a concrete requirement. |
| Plotting | Direct SVG, optional R | Retain existing R scientific plots initially. Simple SVG summaries can follow. |

Upstream references:

- [csv-parser](https://github.com/vincentlaucsb/csv-parser)
- [unordered_dense](https://github.com/martinus/unordered_dense)
- [doctest](https://github.com/doctest/doctest)

Context7 resolution was attempted for all three libraries but returned unrelated
projects. The inspection therefore used the actual vendored source and upstream
documentation. Version numbers above describe the inspected copies, not a claim
that each is the latest release. AntRepCLA's csv.hpp matches its committed SHA-256:
`330b0d8950ff8b566ae726c022e2396f1dcdee36cb9ad6f6c1badc8febbf4319`.

Vendor from identifiable upstream releases/commits. Record upstream URL, version,
commit or release asset, SHA-256, and license for every file. Preserve license
notices and bundled component notices. Builds and tests must not download code.
The selected libraries use MIT licensing; the full vendored notices remain
authoritative. A header reporting a version does not establish release provenance.

## TSV contracts

- Explicit tab delimiter, UTF-8 text, documented header and quoting rules;
  no delimiter inference. Validate required and duplicate header names.
- Set csv-parser's variable-column policy to `THROW`. Its default `IGNORE_ROW`
  can silently discard malformed records. Headerless count inputs also need
  explicit field-count validation.
- Reject missing IDs, duplicate gene/sample identifiers, invalid counts,
  negative/fractional counts, overflow, and unexpected samples. Report filename
  and record information on failure.
- Preserve reference gene IDs as primary keys. Gene symbols are optional display
  annotations; failed or missing annotation must never remove a count row.
- Establish a canonical gene set and fail on unexpected missing genes rather
  than silently applying an inner join. Zero is a measured count, not the default
  interpretation of missing data.
- Preserve declared sample order and deterministic gene order. Never serialize
  hash-table iteration order as a scientific contract.
- Explicitly map sequencing runs to biological samples. Pool technical runs only
  under a declared, validated policy; never infer biological replication from IDs.
- Use a writer with defined escaping, or enforce a restricted output alphabet.
  AntRepCLA's manual tab concatenation is appropriate for its restricted fields,
  but arbitrary sample metadata can contain tabs, quotes, and newlines.
- Parse incrementally where possible. Write completed outputs atomically and
  detect stream errors; a partial file must not count as a completed stage.

## Build and HPC behavior

The native utility build should require a supported `g++` with C++20, GNU Make,
and normal Linux shell tools. It should have `make`, `make test`, and
`make check` targets, generated header dependencies, and explicit parallelism.
This does not remove the external alignment/counting tools, Slurm, reference
assets, or the temporarily retained R environment from the full pipeline.

Do not copy AntRepCLA's unconditional host-wide `nproc` parallelism into an HPC
job. Respect an explicit worker limit and the scheduler allocation. A genuinely
serial parser path should disable parser threading, not merely request one
speculative worker. Avoid nested pools multiplying CPU use.

Make handles file dependencies within an allocation. Slurm arrays and job
dependencies handle distributed scheduling. Account for failed tasks, retries,
configuration/reference changes, and incomplete outputs explicitly. Compile
portable binaries for cluster nodes; do not default to `-march=native`.

## Implementation order and acceptance gates

1. Freeze a validated small count matrix and explicit sample/annotation fixtures.
   Separate corrections to existing behavior from behavior-preserving migration.
2. Add the C++ TSV utility, vendored parser/test headers, and Make build. Accept
   local files first, allowing analysis to run without metadata-service access.
3. Add native count merging and sample validation. Require exact gene IDs, counts,
   sample order, and deterministic outputs; test malformed/quoted/empty fields,
   CRLF, duplicate records, missing genes, integer overflow, and failed writes.
4. Adapt the R input boundary to explicit sample metadata, design, reference,
   and contrast. Keep DESeq2/apeglm pinned and validate sample alignment and outputs.
5. Replace orchestration with Make/Slurm only after resume/failure tests pass.
   Prepare environments and external binaries before compute jobs; remove per-job
   dependency installation.
6. Benchmark alignment/quantification alternatives separately. A Salmon branch
   changes statistical input semantics and is not merely a faster merger.

Integer/count outputs should compare exactly. Statistical outputs need documented
numerical tolerances, matching NA/filter behavior, and stable contrast semantics.
A future native statistical backend additionally needs calibration and power
validation. Keep the reference implementation until those gates pass.

AntRepCLA's shell worker-invariance loop names files by worker count but omits
`--workers` from the invocation. Its unit tests do exercise selected worker counts,
but our end-to-end test must explicitly pass each count being tested.
