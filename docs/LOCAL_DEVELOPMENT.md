# Local development

## Native utilities

```sh
make -j2 check
make -j2 sanitize
build/rnaseq validate-samples tests/fixtures/p0/samples.tsv
build/rnaseq validate-counts tests/fixtures/p0/expected/counts.tsv
build/rnaseq validate-tsv tests/fixtures/p0/runs.tsv
build/rnaseq validate-bundle tests/fixtures/p0/samples.tsv \
  tests/fixtures/p0/runs.tsv tests/fixtures/p0/references.tsv \
  tests/fixtures/p0/analysis.tsv tests/fixtures/p0/contrasts.tsv
build/rnaseq merge-counts tests/fixtures/p0/samples.tsv \
  tests/fixtures/p0/runs.tsv tests/fixtures/p2/genes.tsv \
  tests/fixtures/p2/inputs.tsv /tmp/rnaseq-counts.tsv \
  tests/fixtures/p2/partial-annotation.tsv
build/rnaseq validate-counts-for-samples \
  tests/fixtures/p0/samples.tsv /tmp/rnaseq-counts.tsv
build/rnaseq validate-annotation \
  tests/fixtures/p2/genes.tsv tests/fixtures/p2/partial-annotation.tsv
```

The native build requires GCC 12+ with C++20, GNU Make, Bash and sha256sum;
GCC 15.2.0 is the compiler actually tested here. No R, Python, Slurm, or network
is needed for `make check`. `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, and
`BUILD` can be overridden. Dependency `.d` files track project/vendor headers.
No `-march=native`, OpenMP, automatic host core detection, or parser threads.

The CLI rejects malformed quoting/UTF-8/rows, duplicate or empty headers,
embedded cell separators, invalid sample IDs, duplicate samples/genes, and
invalid/overflowing uint64 counts. Diagnostics identify source, record, and
column. It accepts initial BOM, CRLF, escaped quotes, and empty optional fields.
`validate-tsv` checks table syntax only. Sample/count commands add their stated
schema checks. `validate-bundle` checks all five manifests together: complete
run/sample mapping, library layout/strand consistency, readable regular input
files, reference roles and SHA-256 content, declared covariates, analysis settings,
and categorical contrast levels/references. It resolves each data path against
the directory of the manifest containing it. See [DATA_CONTRACT_V1.md](DATA_CONTRACT_V1.md)
for exact rules. Counts are uint64 at this native boundary; the P3 offline R
interface enforces R's narrower integer limit and model-validity checks.

Reference hashing uses `sha256sum` on PATH through a POSIX subprocess. It passes
the already-opened file through standard input, without a shell or filename
arguments, checks exit status and digest output, and reads files without loading
their contents into native memory. This makes the existing build-time checksum
utility a runtime requirement for bundle validation too. The local implementation
tested is uutils coreutils 0.8.0; its output follows the GNU coreutils interface.
Hashes verify current file contents, not later execution snapshots or R model
validity. The command does not inspect FASTQ/GTF content or run analysis tools.

Small manifests and annotation tables are materialized in memory. P2 count
adapters stream canonical two-column TSV, exact headerless legacy rows, and
featureCounts tables. The input manifest selects formats and the exact
featureCounts count column; the adjacent `.summary` file is mandatory. Merge
uses an explicit reference gene universe, sums technical runs with overflow
checks, and publishes only a fully validated canonical matrix. Its retained
state is O(genes × samples + genes); this is not a constant-memory claim.

Optional annotation is checked separately and never changes counts. See
`tests/fixtures/p2/README.md` for the independently reviewed fixture arithmetic
and provenance. `make check` runs the required P2 oracle comparison. When
retained ignored P0 evidence exists, run
`bash tests/integration/check_p2_retained.sh build/rnaseq` to additionally compare
against the reviewed baseline and real featureCounts 2.0.6 SE/PE outputs. That
optional check does not replace the reproducible fixture gate and does not
establish alignment-backend or statistical parity.

## Analysis dependencies

Preparation uses shell or R with jsonlite; compute uses the native executables
and pinned R/Bioconductor. Python and Snakemake are not required. Mamba is needed
only during explicit online environment creation. For a fresh Linux x86_64 setup:

```sh
mkdir -p .deps tests/output
bash tools/toolchain.sh stage runtime-r --prefix .deps/runtime-r
# The alignment prefix is a staging source, never the compute environment.
bash tools/toolchain.sh stage alignment --prefix .deps/alignment
bash tools/stage_native_runtime.sh .deps/alignment .deps/runtime-tools
.deps/runtime-r/bin/Rscript --vanilla tools/toolchain.R preflight \
  --alignment .deps/runtime-tools --r-prefix .deps/runtime-r \
  --out tests/output/toolchain-preflight
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p4_toolchain.R
```

All destinations must be absent. A failed installation is retained for diagnosis;
compute never installs or fetches dependencies. The shell bootstrap can create R
without a preinstalled R interpreter. For alignment-only bootstrap verification,
use a system R with jsonlite or stage runtime-r first and add its `bin` to PATH.
The R `toolchain.R stage` entry remains available where R/jsonlite already exist.

Preflight checks exact package name/version/build and archive SHA256 identity,
explicit URL/MD5 locks, retained recipes/licenses, the native bundle inventory,
actual versions, pinned R packages, DESeq2/apeglm computation, and PNG output.
Some explicit Mamba installations omit SHA256 in installed records; in that case,
preflight verifies the retained package archive against the pinned SHA256 rather
than accepting MD5 alone. Keep the package cache for this verification.

Preflight constrains auxiliary threads and excludes user R libraries. Successful
output contains `verified.json`, version logs, `R.log`, `MA.png`, and sessionInfo;
failed staging logs remain unpublished for inspection. `--r-kind p0` and an
explicit historical prefix permit historical package checks only. Current defaults
are `.deps/runtime-tools` and `.deps/runtime-r`; see [P7_RUNTIME.md](P7_RUNTIME.md).
On hosts supporting user/network namespaces, prepend `unshare -Urn` for an
independent network-disabled check.

Independent synthetic sequence/count oracles and real tool checks are:

```sh
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p0_fixture.R
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p4_alignment.R
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p0_r.R
```

`check_p0_alignment.R` is retained as an entry-point alias for the broader P4
gate, which checks both native backends, SE/PE, and threads 1/2 against the
independent P0 origins. The fixture generator is `tools/generate_p0_sequences.R`;
`--output DIR` permits regeneration without modifying the repository fixture.
No expected count tables are generated. Historical legacy-bug probes are retired;
recover them only in the isolated tree described in [LEGACY_RECOVERY.md](LEGACY_RECOVERY.md).

Package/source recorders are `tools/record_environment.R PREFIX` and
`tools/record_tool_recipes.R PREFIX [--out DIR]`. They preserve package identities,
archive hashes, original recipes, patches, and license texts. Recipe recording
supports the rendered source blocks of the five pinned tool recipes and rejects
unsupported syntax rather than guessing. Both are preparation tools, not compute.

## P3 offline DESeq2 interface

Use the pinned R environment directly:

```sh
.deps/runtime-r/bin/Rscript --vanilla run_deg_analysis_offline.R \
  --counts tests/fixtures/p3/counts.tsv \
  --samples tests/fixtures/p3/samples.tsv \
  --analysis tests/fixtures/p3/analysis.tsv \
  --contrasts tests/fixtures/p3/contrasts.tsv \
  --annotation tests/fixtures/p3/annotation.tsv \
  --out tests/output/p3/manual \
  --workers 1
```

The output path must be absent and its parent must exist and be writable. Remove
or choose a different manual output path before repeating the command. The script
does not download data, install packages, invoke a shell, or use GEO/BioMart.

Run the complete public-interface suite with:

```sh
mkdir -p tests/output/p3
.deps/runtime-r/bin/Rscript --vanilla tests/integration/check_p3_r.R
```

This performs real DESeq2/apeglm fits, independent-oracle comparisons, workers
1/2/4 parity, failure-preservation checks, and zero/single-significant-gene paths,
so it is slower than the native suite. `P3_SCHEMA_ONLY=1` selects the focused
output-header collision regression but is not the P3 completion gate. Generated
fixture TSVs can be reproduced with
`.deps/runtime-r/bin/Rscript --vanilla tools/generate_p3_fixture.R`; their hashes should
remain unchanged.

Serena's project configuration now requests Bash and R language servers. The R
`languageserver` package version 0.3.19 is installed in the user R library; a
running Serena process may need project reactivation or restart before exposing
the new R symbols.

## Native workflow and real study checks

The local workflow, verified resume, reference cache and explicit Slurm submission
interfaces are documented in [P5_WORKFLOW.md](P5_WORKFLOW.md). `make check` includes
their offline orchestration regressions. `Rscript --vanilla tests/integration/check_p5_real.R`
checks real P0 alignment/merge and the expected insufficient-data R failure.

[P6_REPORT_DATA.md](P6_REPORT_DATA.md) records the TeX-derived GSE80336 inputs,
download checksums, published-count validation, and the user-authorized 30-minute
reduced local check. These opt-in real-data runs are separate from the default
native gate. The 10k-read thinning attempt was too sparse for default normalization;
its retained failure does not justify silently changing the statistical method.
