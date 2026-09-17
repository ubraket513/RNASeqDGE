# P2 count adapters and merge implementation plan

> For agentic workers: execute with the subagent-driven-development or
> executing-plans skill, using test-first changes and an independent review.

Goal: fulfill the approved P2 gate in `IMPLEMENTATION_PLAN.md` without changing
the active HISAT2 workflow or R statistical behavior. Work in the current `main`
checkout as explicitly authorized; preserve uncommitted P1 work.

Spec: `IMPLEMENTATION_PLAN.md` (counts/annotation and P2) plus
`DATA_CONTRACT_V1.md`. The following interface decisions make P2 concrete.

## Interfaces and decisions

- `merge-counts SAMPLES RUNS GENES INPUTS OUTPUT [ANNOTATION]` writes a canonical
  integer matrix. `GENES` is a one-column `gene_id` reference universe, unique,
  nonempty, sorted bytewise on output. Its provenance must tie it to the reference;
  P2 does not infer the universe from the first run or gene symbols.
- `INPUTS` has exactly `run_id`, `format`, `path`, `column`. It covers every run
  exactly once, no extras. Mapping to samples comes only from `RUNS`.
- Formats: `tsv` requires `gene_id,count` header; `legacy` is exactly two
  headerless columns; `featurecounts` requires the six standard annotation
  columns plus count columns. `column` is empty for tsv/legacy and names an
  exact count header for featurecounts. No automatic format/column inference.
- FeatureCounts inputs require the accompanying `PATH.summary`, with unique
  nonempty status rows, a selected column matching the count table and an
  `Assigned` row. All counts are unsigned integers. Preserve both source files.
  Summary assignments are not assumed equal to summed gene counts because
  upstream counting policies can permit multiple assignments.
- Every run must contain exactly the reference universe, each gene once.
  Keep measured zero rows, original IDs/suffixes, and sample manifest order.
  Add technical-run counts with overflow checks; never fill a missing gene/sample
  with zero. `validate-counts-for-samples SAMPLES COUNTS` checks exact columns.
- Optional annotation is independently validated against the reference universe:
  exact `gene_id,gene_symbol,biotype,chromosome` schema, unique gene IDs, subset
  allowed, empty optional values allowed, duplicate symbols allowed. It is not
  joined into counts. `validate-annotation GENES ANNOTATION` exposes that check.
- Input paths use their containing manifest directory, preserving symlink/`..`
  semantics and allowing absolute paths. Merge checks mapping/library metadata;
  archived FASTQ files need not remain present once counts exist.
- Stream count input rows. Retain the gene index and numeric aggregate matrix:
  memory is O(genes × samples + genes), not constant-memory matrix processing.
  No whole count input text/cell tables, worker threads, or new dependencies.
- Validate all inputs before publishing. Write to an exclusive temporary file
  alongside OUTPUT, check write/close errors, then atomically rename. On failure
  preserve prior output, remove owned temporary files, reject output aliases to
  any input or summary. No multiprocess resume/cache/locking claim (P5).

## Task 1: shared streaming TSV records

Files: `include/rnaseq/table.hpp`, `src/table.cpp`,
`tests/unit/test_table.cpp`.

Public interface: `TsvReader(std::istream&, std::string source,
bool leading_comments=false)`, `bool next(std::vector<std::string>&)`,
`std::size_t record() const`. It reads one physical record at a time, keeps the
existing strict UTF-8/quoting/BOM/CRLF/width policy, and optionally skips only
leading featureCounts `#` lines. `read_table` builds small tables from this reader.

- [ ] Add runtime regression tests for lazy row consumption, late malformed rows,
  CRLF/quoted/UTF-8 cells, headerless rows, and leading-only comments.
- [ ] Verify failures with an unimplemented reader, then implement and pass all
  existing table tests plus new tests.

## Task 2: adapters and checked aggregation

Files: `include/rnaseq/counts.hpp`, `src/counts.cpp`,
`tests/unit/test_counts.cpp`. Consume Task 1 reader and existing count helpers.
Expose `merge_counts(samples,runs,genes,inputs,output,annotation={})`,
`validate_counts_for_samples(samples,counts)`,
`validate_annotation(genes,annotation)` with filesystem paths.

- [ ] Add tests with hand-calculated two-sample/multi-run values, shuffled genes,
  all formats, explicit featureCounts column selection and required summaries.
- [ ] Verify failure before implementing merge, then implement. Test absent and
  duplicate runs/genes, universe mismatch, invalid/fractional/overflow counts,
  annotation independence, path semantics, and output preservation/alias checks.
- [ ] Reuse P1 run validation by factoring its mapping checks into
  `validate_run_mapping(const Table&, const Table&, bool check_fastq_files=true)`;
  P1 retains full path checks, P2 passes false but still validates required path
  cells/layout/strand and complete sample mapping. Parent owns this factoring.

## Task 3: CLI, reproducible fixtures, review and gate

Files: `src/main.cpp`, `Makefile`, `tests/integration/check_counts.sh`,
`tests/fixtures/p2/`, `docs/{DATA_CONTRACT_V1,LOCAL_DEVELOPMENT,
IMPLEMENTATION_PROGRESS,SESSION_HANDOFF}.md`. Parent owns these files.

- [ ] Add CLI dispatch/build entries and an offline shell integration gate.
- [ ] Compare P0 count inputs to its unchanged hand-derived oracle byte-for-byte;
  compare the generated P0 reviewed baseline when present, clearly distinguish
  optional historical evidence from required reproducible tests.
- [ ] Test featureCounts inputs/summary with independently documented synthetic
  fixtures and the real retained P0 alignment outputs when present.
- [ ] Independent review, fix material findings with regression tests.
- [ ] Run `make -j2 check`, `make -j2 sanitize`, `git diff --check`; record actual
  results and remaining P3/P4/P5/P6 limits. Do not claim R parity or performance.

## Preflight interface review

| Tasks | Shared surface | Resolution |
| --- | --- | --- |
| 1 / 2 | TsvReader rows and physical record numbers | fixed API above; Task 2 can compile after Task 1 lands |
| 2 / 3 | three path-based public functions and CLI order | fixed API above; parent integrates after declarations |
| P1 / 2 | run validation | preserve P1 behavior through default true; merge permits archived FASTQs |
| 1 | strict existing TSV semantics vs streaming | keep all old tests; no format guessing |
| 2 | arbitrary input order vs sorted output | explicit universe plus checked numeric matrix |
| 3 | required offline gate vs ignored historical outputs | hand oracle mandatory; retained P0 evidence checked separately |
