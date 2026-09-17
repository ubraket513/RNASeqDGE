# Local development

## Native utilities

```sh
make -j2 check
make -j2 sanitize
build/rnaseq validate-samples tests/fixtures/p0/samples.tsv
build/rnaseq validate-counts tests/fixtures/p0/expected/counts.tsv
build/rnaseq validate-tsv tests/fixtures/p0/runs.tsv
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
schema checks. This is not yet a full analysis-bundle validator: run mapping,
reference hash matching, covariate declarations, and model checks remain later
work. Counts are uint64 here; R's narrower integer limit must be enforced at P3.

P1 currently materializes input tables in memory. This is appropriate for small
manifests and fixtures; P2 must add streaming count processing before claiming
large-matrix memory performance. The CLI validates but never rewrites inputs.

## Analysis dependencies

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
