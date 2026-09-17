#!/usr/bin/env bash
set -euo pipefail
binary=${1:?native binary required}
"$binary" validate-samples tests/fixtures/p0/samples.tsv
"$binary" validate-counts tests/fixtures/p0/expected/counts.tsv
for name in runs references contrasts analysis annotation; do
    "$binary" validate-tsv "tests/fixtures/p0/$name.tsv"
done
"$binary" validate-bundle \
    tests/fixtures/p0/samples.tsv \
    tests/fixtures/p0/runs.tsv \
    tests/fixtures/p0/references.tsv \
    tests/fixtures/p0/analysis.tsv \
    tests/fixtures/p0/contrasts.tsv
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT
if "$binary" > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: missing command accepted' >&2
    exit 1
else
    status=$?
fi
test "$status" -eq 2
printf 'gene_id\ts\ngene\t18446744073709551616\n' > "$temporary/overflow.tsv"
if "$binary" validate-counts "$temporary/overflow.tsv" > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: overflow accepted' >&2
    exit 1
fi
grep -q 'record 2, column s' "$temporary/error"
if "$binary" validate-tsv "$temporary/missing.tsv" > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: missing file accepted' >&2
    exit 1
fi
if "$binary" unknown tests/fixtures/p0/samples.tsv > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: unknown command accepted' >&2
    exit 1
fi
printf 'key\tvalue\nversion\t1\ndesign_terms\tcondition\nalpha\t0.05\nfilter\tzero_total\nshrinkage\tapeglm\nunknown\tvalue\n' > "$temporary/analysis.tsv"
if "$binary" validate-bundle \
    tests/fixtures/p0/samples.tsv \
    tests/fixtures/p0/runs.tsv \
    tests/fixtures/p0/references.tsv \
    "$temporary/analysis.tsv" \
    tests/fixtures/p0/contrasts.tsv > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: unknown analysis setting accepted' >&2
    exit 1
fi
grep -q 'analysis.tsv: record 7, column key' "$temporary/error"
echo 'PASS: native CLI checks'
if [ -e /dev/full ]; then
    if "$binary" validate-samples tests/fixtures/p0/samples.tsv > /dev/full 2> "$temporary/error"; then
        echo 'FAIL: stdout write failure accepted' >&2
        exit 1
    fi
fi
