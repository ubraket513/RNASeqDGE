#!/usr/bin/env bash
set -euo pipefail
binary=$(realpath "${1:?native binary required}")
fixture=$(realpath tests/fixtures)
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT

"$binary" merge-counts "$fixture/p0/samples.tsv" "$fixture/p0/runs.tsv" \
    "$fixture/p2/genes.tsv" "$fixture/p2/inputs.tsv" "$temporary/counts.tsv"
cmp "$fixture/p0/expected/counts.tsv" "$temporary/counts.tsv"
"$binary" validate-counts-for-samples "$fixture/p0/samples.tsv" "$temporary/counts.tsv"
"$binary" validate-annotation "$fixture/p2/genes.tsv" "$fixture/p2/partial-annotation.tsv"

# Run from an unrelated directory: paths are relative to their own manifest.
(
    cd "$temporary"
    "$binary" merge-counts "$fixture/p0/samples.tsv" "$fixture/p0/runs.tsv" \
        "$fixture/p2/genes.tsv" "$fixture/p2/mixed-inputs.tsv" mixed.tsv \
        "$fixture/p2/partial-annotation.tsv"
)
cmp "$fixture/p0/expected/counts.tsv" "$temporary/mixed.tsv"

# Failures preserve the prior output and do not modify raw counts or summaries.
cp -R "$fixture/p0" "$temporary/p0"
cp -R "$fixture/p2" "$temporary/p2"
printf 'gene_id\tcount\ngene_A.1\t1\ngene_B.2\t2\ngene_C\t3\n' \
    > "$temporary/p0/counts/run_c2.tsv"
printf 'previous complete output\n' > "$temporary/protected.tsv"
cp "$temporary/protected.tsv" "$temporary/prior.tsv"
sha256sum "$temporary/p2/counts/"* > "$temporary/hashes.before"
if "$binary" merge-counts "$temporary/p0/samples.tsv" "$temporary/p0/runs.tsv" \
    "$temporary/p2/genes.tsv" "$temporary/p2/mixed-inputs.tsv" \
    "$temporary/protected.tsv" > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: missing gene accepted' >&2
    exit 1
fi
cmp "$temporary/prior.tsv" "$temporary/protected.tsv"
sha256sum "$temporary/p2/counts/"* > "$temporary/hashes.after"
cmp "$temporary/hashes.before" "$temporary/hashes.after"

printf 'gene_id\tcontrol_1\ttreated_2\ttreated_1\tcontrol_2\ng\t1\t2\t3\t4\n' \
    > "$temporary/reordered.tsv"
if "$binary" validate-counts-for-samples "$fixture/p0/samples.tsv" \
    "$temporary/reordered.tsv" > "$temporary/out" 2> "$temporary/error"; then
    echo 'FAIL: reordered sample columns accepted' >&2
    exit 1
fi

for command in merge-counts validate-counts-for-samples validate-annotation; do
    if "$binary" "$command" > "$temporary/out" 2> "$temporary/error"; then
        echo "FAIL: missing arguments accepted for $command" >&2
        exit 1
    else
        status=$?
    fi
    test "$status" -eq 2
done
echo 'PASS: count adapters match hand-derived P0 matrix; failure preserves output'
