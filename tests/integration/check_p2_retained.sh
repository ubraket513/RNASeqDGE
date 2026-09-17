#!/usr/bin/env bash
set -euo pipefail
binary=$(realpath "${1:-build/rnaseq}")
fixture=$(realpath tests/fixtures)
alignment=$(realpath tests/output/p0/alignment)
reviewed=$(realpath tests/output/p0/reviewed/counts.tsv)
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT
test -f "$reviewed"

"$binary" merge-counts "$fixture/p0/samples.tsv" "$fixture/p0/runs.tsv" \
    "$fixture/p2/genes.tsv" "$fixture/p2/mixed-inputs.tsv" "$temporary/counts.tsv"
cmp "$reviewed" "$temporary/counts.tsv"
cmp "$fixture/p0/expected/counts.tsv" "$temporary/counts.tsv"

printf 'sample_id\tcondition\nPE\tsynthetic\nSE\tsynthetic\n' > "$temporary/samples.tsv"
printf 'run_id\tsample_id\tfastq_1\tfastq_2\tlayout\tstrandedness\npaired\tPE\tarchived_1.fastq\tarchived_2.fastq\tpaired\tunstranded\nsingle\tSE\tarchived.fastq\t\tsingle\tunstranded\n' \
    > "$temporary/runs.tsv"
printf 'run_id\tformat\tpath\tcolumn\n' > "$temporary/inputs.tsv"
for run in single paired; do
    counts="$alignment/$run.counts.txt"
    test -f "$counts.summary"
    column=$(awk -F '\t' '$1 == "Geneid" {print $7; exit}' "$counts")
    test -n "$column"
    printf '%s\tfeaturecounts\t%s\t%s\n' "$run" "$counts" "$column" >> "$temporary/inputs.tsv"
done
sha256sum "$alignment/"*.counts.txt "$alignment/"*.counts.txt.summary > "$temporary/before"
"$binary" merge-counts "$temporary/samples.tsv" "$temporary/runs.tsv" \
    "$fixture/p2/genes.tsv" "$temporary/inputs.tsv" "$temporary/aligned.tsv"
# Independently reviewed P0 read-origin oracle, not adapter-generated expectations.
printf 'gene_id\tPE\tSE\ngene_A.1\t1\t2\ngene_B.2\t0\t1\ngene_C\t0\t0\ngene_zero\t0\t0\n' \
    > "$temporary/expected.tsv"
cmp "$temporary/expected.tsv" "$temporary/aligned.tsv"
sha256sum "$alignment/"*.counts.txt "$alignment/"*.counts.txt.summary > "$temporary/after"
cmp "$temporary/before" "$temporary/after"
echo 'PASS: reviewed P0 baseline and retained real SE/PE featureCounts outputs'
