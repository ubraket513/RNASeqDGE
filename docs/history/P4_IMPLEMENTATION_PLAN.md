# P4: pinned tools and offline alignment backends

Execute the P4 gate in IMPLEMENTATION_PLAN.md on the user-authorized existing
main checkout. Preserve P1–P3 work. No commits, cluster submissions, or backend
default change. HISAT2 remains baseline; STAR is an explicit option.

## Task 1: native backend runner

Implement C++20 argv-based subprocess execution and explicit offline index and
alignment/counting commands. Reuse existing manifest/validation/hash machinery
where appropriate. Expose a documented CLI with explicit executable directory,
reference FASTA/GTF, reads/layout/strandedness, backend, threads and absent output
directory. Default worker budget is one. Reject invalid budgets and allocations
above SLURM_CPUS_PER_TASK. Invoke tools sequentially so aligner and samtools sort
never each consume a concurrent whole budget; samtools auxiliary threads must
account for its main thread. Bound sort memory explicitly. No shell interpolation,
host core detection, downloads, installation, or arbitrary extra tool arguments.

Use installed pins: HISAT2 2.2.1, STAR 2.7.10b, samtools 1.17 and featureCounts
2.0.6. Keep wrappers as supplied (HISAT2 needs Python/Perl). STAR produces unsorted
BAM followed by samtools sort; preserve BAM, alignment logs, featureCounts raw
output and summary. Explicit paired fragment flags and strandedness mapping;
MAPQ 0, exon/gene_id counting, default exclusion of ambiguous/multimapping reads.
Index policy must be explicit and recorded, with configurable tiny-genome STAR
index parameters rather than inappropriate full-human defaults for tiny fixtures.
Validate all inputs before launch, stage sibling output and publish on success;
preserve logs on failure without publishing successful completion. Propagate
tool exit failures and interruption, terminate/reap active children. P5 owns cache,
run orchestration, resume and cross-run locking.

Own src/process.cpp, include/rnaseq/process.hpp, src/alignment.cpp,
include/rnaseq/alignment.hpp, src/main.cpp, Makefile, tests/unit/test_process.cpp,
tests/integration/check_p4_native.sh, tests/integration/check_p4_alignment.py,
and tests/fixtures/p4 if necessary. Keep native make check offline and independent
of Python/R/alignment tool installations. Add meaningful process/failure/resource
tests with stub executables plus separate real STAR/HISAT2 SE/PE tests using P0
fixtures and independently stated counts/junction/summary expectations. Check
threads 1 and 2 where supported. Verify make check, sanitizer on new native
boundary, and real tool gate. Record any scientific mismatch rather than changing
the oracle to fit output. Update docs/P4_BACKENDS.md with CLI and test evidence.

## Task 2: tool provenance and staging

Controller prepares source/build provenance for pinned executable packages and
transitive dependencies from installed conda records/recipes, an explicit staging
entry point using existing locks, and an offline preflight that checks tool
versions plus pinned R packages, a tiny DESeq2/apeglm computation and PNG device.
Clearly distinguish binary conda staging from source builds and offline compute
from offline redistribution. Retain source URL/hash/license/build flags from
recipes without inventing missing metadata. Never install during preflight.

## Completion

Review new P4 changes, fix concrete findings with focused regressions, update
progress/handoff/development docs and Serena memories. Reuse successful test
evidence for unchanged code; do not rerun P3 statistical parity unless P3 changes.
