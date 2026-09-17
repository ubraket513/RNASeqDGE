#!/usr/bin/env bash
set -euo pipefail
exe=$(realpath "${1:-build/rnaseq}")
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/bin"
cat > "$tmp/bin/hisat2-build-s" <<'STUB'
#!/bin/sh
if [ "$1" = --version ]; then echo 'hisat2-build-s version 2.2.1'; exit; fi
printf '%s\n' "$@" > "$P4_ARGS"
exit 17
STUB
chmod +x "$tmp/bin/hisat2-build-s"
export P4_ARGS="$tmp/args"
base=(index --backend hisat2 --bin-dir "$tmp/bin" --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf --output "$tmp/output")
for budget in 0 -1 two 2147483648; do
  if "$exe" "${base[@]}" --threads "$budget" >"$tmp/log" 2>&1; then exit 1; fi
  test ! -e "$P4_ARGS"
done
if SLURM_CPUS_PER_TASK=1 "$exe" "${base[@]}" --threads 2 >"$tmp/log" 2>&1; then exit 1; fi
test ! -e "$P4_ARGS"
set +e
"$exe" "${base[@]}" >"$tmp/log" 2>&1
status=$?
set -e
test "$status" = 17
test ! -e "$tmp/output"
test -f "$P4_ARGS"
test "$(sed -n 3p "$P4_ARGS")" = -p
test "$(sed -n 4p "$P4_ARGS")" = 1
grep -q '17' "$tmp/log"
test "$(find "$tmp" -name index.log | wc -l)" = 1
printf 'PASS P4 native: resource validation, tool exit propagation, failure staging\n'
# A real child process must be killed and reaped on interruption, even if it
# ignores TERM. The PID marker is a synchronization boundary, not a timing guess.
cat > "$tmp/bin/hisat2-build-s" <<'STUB'
#!/bin/sh
if [ "$1" = --version ]; then echo 'hisat2-build-s version 2.2.1'; exit; fi
trap '' TERM INT
printf '%s\n' "$$" > "$P4_ARGS"
while :; do sleep 1; done
STUB
rm "$P4_ARGS"
"$exe" "${base[@]}" >"$tmp/interrupt.log" 2>&1 &
runner=$!
for _ in $(seq 1 200); do test -f "$P4_ARGS" && break; sleep .01; done
test -f "$P4_ARGS"
child=$(cat "$P4_ARGS")
kill -TERM "$runner"
set +e
wait "$runner"
status=$?
set -e
test "$status" = 143
! kill -0 "$child" 2>/dev/null
test ! -e "$tmp/output"
test "$(find "$tmp" -name COMPLETE | wc -l)" = 0
printf 'PASS P4 native: interruption propagates and active child is reaped\n'

# A successful stub index exercises literal paths and deterministic parameters.
cat > "$tmp/bin/hisat2-build-s" <<'STUB'
#!/bin/sh
if [ "$1" = --version ]; then echo 'hisat2-build-s version 2.2.1'; exit; fi
printf '%s\n' "$@" > "$P4_ARGS"
for prefix; do :; done
for i in 1 2 3 4 5 6 7 8; do printf x > "$prefix.$i.ht2"; done
STUB
literal="$tmp/space ; \$(touch NEVER)"
"$exe" index --backend hisat2 --bin-dir "$tmp/bin" --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf --output "$literal" --threads 2
test -f "$literal/COMPLETE"
test "$(sed -n 4p "$P4_ARGS")" = 2
test ! -e NEVER
# Preflight failures may not launch even a version command.
rm "$P4_ARGS"
if "$exe" index --backend hisat2 --bin-dir "$tmp/bin" --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf --output "$literal" >"$tmp/existing.log" 2>&1; then exit 1; fi
test ! -e "$P4_ARGS"
printf 'PASS P4 native: explicit budget, literal paths, staged publication, existing-output rejection\n'
# A competing publisher must not result in a successful-completion marker in
# the failed stage. A stub simulates creation of output during tool execution.
cat >> "$tmp/bin/hisat2-build-s" <<'STUB'
mkdir "$P4_COMPETING_OUTPUT"
STUB
export P4_COMPETING_OUTPUT="$tmp/competing"
if "$exe" index --backend hisat2 --bin-dir "$tmp/bin" --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf --output "$P4_COMPETING_OUTPUT" >"$tmp/competing.log" 2>&1; then exit 1; fi
test "$(find "$tmp" -path '*competing.staging-*' -name COMPLETE | wc -l)" = 0
# Exact pins must reject version prefixes before any index work starts.
cat > "$tmp/bin/hisat2-build-s" <<'STUB'
#!/bin/sh
if [ "$1" = --version ]; then echo "hisat2-build-s version $P4_STUB_VERSION"; exit; fi
printf work > "$P4_ARGS"
for prefix; do :; done
for i in 1 2 3 4 5 6 7 8; do printf x > "$prefix.$i.ht2"; done
STUB
rm -f "$P4_ARGS"
for version in 2.2.10 12.2.1 2.2.1-suffix; do
  if P4_STUB_VERSION="$version" "$exe" index --backend hisat2 --bin-dir "$tmp/bin" --fasta tests/fixtures/p0/reference.fa --gtf tests/fixtures/p0/reference.gtf --output "$tmp/near-$version" >"$tmp/near.log" 2>&1; then
    echo "FAIL: accepted near-match version $version" >&2; exit 1
  fi
  test ! -e "$P4_ARGS"
  grep -q 'expected pinned version 2.2.1' "$tmp/near.log"
done
printf 'PASS P4 native: exact version pins reject near matches before work\n'
