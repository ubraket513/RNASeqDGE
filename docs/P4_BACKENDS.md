# P4 offline alignment backends

The native C++20 CLI builds a genome-only index and runs one alignment/counting
job. Supply an installed executable directory; execution never installs or
fetches tools. Required pins are HISAT2 2.2.1, STAR 2.7.10b, samtools 1.17,
and featureCounts 2.0.6. Version output is checked and retained. HISAT2 launchers
remain upstream wrappers; their Python/Perl runtime dependencies must be present.

From the repository root, for the P0 tiny genome:

```sh
make
build/rnaseq index --backend star --bin-dir .deps/alignment/bin \
  --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf \
  --star-sa-bases 3 --star-chr-bits 10 --output /tmp/p4-star-index
build/rnaseq align-count --backend star --bin-dir .deps/alignment/bin \
  --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf \
  --index /tmp/p4-star-index --reads1 tests/fixtures/p0/reads/single.fastq \
  --layout single --strandedness unstranded --output /tmp/p4-star-single
```

For HISAT2 use `--backend hisat2` in both commands and omit the two STAR index
options. For PE use `--layout paired --reads1 R1.fastq --reads2 R2.fastq`.
`--strandedness` is required: `unstranded`, `forward`, or `reverse` maps to
featureCounts `-s 0`, `-s 1`, or `-s 2`. HISAT2 additionally uses F/R for SE and
FR/RF for PE. STAR alignment stays strand-agnostic; counting applies the specified
library strand. P4 accepts readable uncompressed FASTQ files only; `.gz` inputs
are rejected. This is file-level validation, not a FASTQ sequence integrity check.

`--threads N` defaults to **1**. It must be a positive integer within
`SLURM_CPUS_PER_TASK` when that variable is set. The CLI never detects host cores.
Each subprocess completes before the next starts. `samtools sort` receives
`-@ N-1 -m 256M`, reserving its main thread within the worker budget and explicitly
bounding approximate memory per sorting thread. This is not a total process RSS
limit; alignment/index memory remains tool-dependent.

The index policy is explicitly **genome-only**, matching P0 HISAT2 discovery.
GTF is not injected into either index; it supplies exon annotations during
counting. STAR defaults are SA bases 14 and chromosome bits 18; choose appropriate
values explicitly for small references (P0 uses 3 and 10). Allowed ranges are
1–14 and 1–18 respectively. Index metadata records those values, backend,
FASTA/GTF SHA-256 hashes, and worker budget. Alignment rejects a mismatched
reference/backend identity. HISAT2's small `.ht2` index format is supported in P4;
large `.ht2l` indexes are outside this gate.

STAR writes unsorted BAM, followed by sequential samtools sorting; HISAT2 writes
SAM, followed by the same sorting boundary. Both retain the sorted `aligned.bam`,
original alignment, tool logs, version logs, argv records (`commands.tsv`), raw
`counts.txt`, and `counts.txt.summary`. Counting is MAPQ 0, exon/gene_id, with
ambiguous and multimapping reads excluded by featureCounts defaults. PE explicitly
uses both `-p` and `--countReadPairs`, so one fragment contributes one count.
No arbitrary extra tool arguments are accepted.

The output directory must not exist and its parent must exist. Work runs in a
unique sibling `.staging-*` directory, then publishes by rename on success.
Failures preserve that stage and logs and report its path; tool exit status is
propagated. SIGINT/SIGTERM kill the active subprocess group, reap the direct child,
and return 130/143. `COMPLETE` marks success. This boundary does not implement
cache lookup, resume, cross-run locks, or orchestration; these belong to P5.
Recorded argv and featureCounts BAM column names reflect the execution-stage
paths; those historical paths remain useful provenance after publication.

## Verification

`make check` is offline and needs neither Python/R nor alignment installations.
It covers argv literals, failure status, missing executables, resource rejection,
default worker budget, literal output paths, publication, and interruption/reaping
through small shell stubs. Existing P1/P2 checks remain in this gate.

The separately installed real-tool gate is:

```sh
python3 tests/integration/check_p4_alignment.py
```

It builds both indexes at one and two threads, and checks SE/PE jobs for each.
The independent P0 origins require SE counts A=2, B=1, C=0, zero=0, three assigned
reads and one ambiguous read; PE counts A=1 and all others zero, with one assigned
fragment. SE junction CIGAR must be `37M200N38M`. These are scientific expectations,
not captured native output. All eight combinations passed on 2026-09-17; evidence
is in `tests/output/p4-real-smg3k40j/verified.json` with complete outputs/logs.

Context7 samtools documentation was consulted for sorting. The implementation
uses the installed pinned CLI options, including samtools auxiliary threads and
featureCounts paired-fragment semantics.

The final real gate also passed inside a network-disabled user namespace:
`unshare -Urn python3 tests/integration/check_p4_alignment.py`, with transcript
`tests/output/p4-network-isolated.log` and evidence
`tests/output/p4-real-hhh9ogyf/verified.json`. The full native check passed with
AddressSanitizer/UBSan using `-fno-pie`/`-no-pie`; after the final publication fix,
the rebuilt sanitized CLI passed the focused native P4 gate again. Publication
uses Linux `renameat2(RENAME_NOREPLACE)` and removes the stage's completion marker
if a competitor creates the output before publication. P4 therefore requires
Linux with `renameat2` support; this check is distinct from P5 cross-run locking.

The reviewed subprocess boundary prepends `--bin-dir` to each child's PATH, so
upstream `/usr/bin/env python` and Perl wrappers resolve staged interpreters first.
It builds a private child environment; the parent PATH is unchanged. Paths with
spaces or shell metacharacters stay literal; colons in the executable directory
are rejected because PATH uses colons as separators. Version checks require exact
tool-specific version lines, rejecting near versions such as HISAT2 2.2.10.
After those fixes, all eight cases passed again with networking disabled and
ambient `PATH=/usr/bin:/bin`: see `tests/output/p4-network-isolated-review.log`
and `tests/output/p4-real-muu30lsa/verified.json`. Native checks passed 36 cases/
222 assertions; focused sanitized process tests and the P4 native gate also passed.
