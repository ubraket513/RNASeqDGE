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

P4 adds an explicit staging/preflight boundary (Python 3.11+ standard library;
Mamba is required only for staging). To prepare a new Linux x86_64 prefix:

```sh
python3 tools/toolchain.py stage alignment --prefix .deps/alignment-new
python3 tools/toolchain.py stage p0 --prefix .deps/r-new
```

These commands install the exact existing explicit locks and verify installed
name/version/build/archive-SHA-256 records. Prefixes must be absent. Staging is
online; a failed install is retained for diagnosis and must not be treated as
ready. Existing prefixes are never overwritten. Compute never calls staging.

For an installed environment, run the offline preflight:

```sh
mkdir -p tests/output
python3 tools/toolchain.py preflight --out tests/output/toolchain-preflight
python3 tests/integration/check_p4_toolchain.py
```

The output directory must be absent. `--alignment` and `--r-prefix` select other
prefixes. Preflight checks both full package locks, retained recipe/license hashes,
actual STAR/HISAT2/build/samtools/featureCounts versions, pinned R packages,
DESeq2/apeglm computation, and a PNG device. It constrains auxiliary threads and
excludes user R libraries. Successful output includes `verified.json`, version
logs, `R.log`, `MA.png` and `sessionInfo.txt`; failed temporary output is retained
without publishing the requested directory. Package metadata checks do not prove
every installed file is unmodified; actual probed entry-point hashes are recorded.

On Linux hosts allowing user/network namespaces, prefix the command with
`unshare -Urn` to independently verify execution without networking. This passed
locally for both the existing and a newly reconstructed alignment prefix.
See [vendor toolchain provenance](../vendor/toolchain/README.md) for source/build
records and the distinction between locked binary deployment and source rebuilds.
Native backend CLI and real-tool checks are in [P4_BACKENDS.md](P4_BACKENDS.md).

Installed locally under ignored `.deps/`:

- `.deps/p0`: Python/Snakemake plus R/Bioconductor and required plotting packages.
- `.deps/alignment`: legacy HISAT2 2.2.1, Subread 2.0.6, samtools 1.17, STAR 2.7.10b.

Separate prefixes avoid the old alignment stack's zlib conflict with current R.
Explicit lock files in `config/` pin all resolved package URLs, builds and MD5s;
the accompanying provenance JSON records SHA-256, licenses, dependencies, and
post-link R source URLs/MD5s from the installed recipes. The YAML files describe requested constraints and
are **not** reproducibility locks. Use the explicit files to reconstruct:

```sh
MAMBA_ROOT_PREFIX="$PWD/.deps/mamba" CONDA_PKGS_DIRS="$PWD/.deps/pkgs" \
  mamba create --prefix "$PWD/.deps/p0" --file config/p0-linux-64.explicit.txt --yes
MAMBA_ROOT_PREFIX="$PWD/.deps/mamba" CONDA_PKGS_DIRS="$PWD/.deps/pkgs-align" \
  mamba create --prefix "$PWD/.deps/alignment" --file config/alignment-linux-64.explicit.txt --yes
```

Environment creation is an online staging step. Bioconductor data packages can
download source archives in post-link scripts. The org.Hs.eg.db installer verified
its recorded MD5 and removed the tarball afterward; no archive SHA-256 was
captured for that extra download. This is not an offline installation bundle.
Compute/test execution performs no package installation or metadata downloads.
Do not combine these prefixes through `LD_LIBRARY_PATH`; invoke installed
executables using the PATH below. A scheduler is not installed on this machine.

```sh
export PATH="$PWD/.deps/p0/bin:$PWD/.deps/alignment/bin:$PATH"
export XDG_CACHE_HOME="$PWD/.deps/cache"
python tests/integration/check_p0_fixture.py
python tests/integration/check_p0_legacy.py
python tests/integration/check_p0_alignment.py
Rscript --vanilla tests/integration/check_p0_r.R
```

All use synthetic local inputs. The legacy tests assert observed bugs separately
from the reviewed integer oracle; the R test is backend functionality, not
original-study statistical parity. R's parametric dispersion fit may fall back
to a local regression on this synthetic dataset; the log records that choice.
R outputs include sessionInfo, raw/shrunken tables, and a real PNG.

Raw outputs/logs live in ignored `tests/output/p0/`. No production config, aligner
default, or submission script is changed by these checks. The two existing
cluster smoke-test scripts still need their cluster tools and online study data;
they are not replaced by these local checks. SRA Toolkit and Slurm were not
installed because the completed checks use local reads and no cluster jobs.

## P3 offline DESeq2 interface

Use the pinned R environment directly:

```sh
.deps/p0/bin/Rscript --vanilla run_deg_analysis_offline.R \
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
.deps/p0/bin/Rscript --vanilla tests/integration/check_p3_r.R
```

This performs real DESeq2/apeglm fits, independent-oracle comparisons, workers
1/2/4 parity, failure-preservation checks, and zero/single-significant-gene paths,
so it is slower than the native suite. `P3_SCHEMA_ONLY=1` selects the focused
output-header collision regression but is not the P3 completion gate. Generated
fixture TSVs can be reproduced with
`.deps/p0/bin/Rscript --vanilla tools/generate_p3_fixture.R`; their hashes should
remain unchanged.

Serena's project configuration now requests Bash and R language servers. The R
`languageserver` package version 0.3.19 is installed in the user R library; a
running Serena process may need project reactivation or restart before exposing
the new R symbols.

## Native workflow and real study checks

The local workflow, verified resume, reference cache and explicit Slurm submission
interfaces are documented in [P5_WORKFLOW.md](P5_WORKFLOW.md). `make check` includes
their offline orchestration regressions. `python3 tests/integration/check_p5_real.py`
checks real P0 alignment/merge and the expected insufficient-data R failure.

[P6_REPORT_DATA.md](P6_REPORT_DATA.md) records the TeX-derived GSE80336 inputs,
download checksums, published-count validation, and the user-authorized 30-minute
reduced local check. These opt-in real-data runs are separate from the default
native gate. The 10k-read thinning attempt was too sparse for default normalization;
its retained failure does not justify silently changing the statistical method.
