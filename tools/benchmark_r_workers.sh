#!/usr/bin/env bash
# Sequential, memory-guarded full-study benchmark. Requires Linux procfs and GNU time.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
input=$(realpath "${1:-$repo/tests/output/p6-inputs/full36}")
output=${2:-$repo/tests/output/p7-benchmark-$(date -u +%Y%m%dT%H%M%SZ)}
rscript=$(realpath "${3:-$repo/.deps/p0/bin/Rscript}")
baseline=${4:-$repo/tests/output/p6-full36}
mkdir "$output"
output=$(realpath "$output")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset R_HOME R_LIBS LD_LIBRARY_PATH LD_PRELOAD
export R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null
reserve_kb=${P7_RESERVE_KB:-655360}
limit_seconds=${P7_BENCHMARK_SECONDS:-1800}
start=$(date +%s)
active=
cleanup() { if [ -n "$active" ]; then kill -TERM -- "-$active" 2>/dev/null || true; sleep 1; kill -KILL -- "-$active" 2>/dev/null || true; fi; }
trap cleanup EXIT
trap 'cleanup; exit 130' INT
trap 'cleanup; exit 143' TERM
printf 'workers\tstatus\twall_seconds\tpeak_group_rss_kb\tmin_available_kb\tinitial_available_kb\trequired_headroom_kb\n' > "$output/benchmark.tsv"
printf 'epoch\tworkers\tgroup_rss_kb\tavailable_kb\n' > "$output/memory.tsv"
printf 'reserve_kb\t%s\ntime_limit_seconds\t%s\n' "$reserve_kb" "$limit_seconds" > "$output/guard.txt"
one_peak=0
for workers in 1 2 4; do
 available=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
 # An observed single-worker peak plus conservative extra copy allowance.
 required=$((reserve_kb + 1100000))
 if [ "$one_peak" -gt 0 ]; then required=$((reserve_kb + one_peak + (workers - 1) * one_peak / 2)); fi
 if [ "$available" -lt "$required" ]; then
  printf '%s\tinsufficient_headroom\t0\t0\t%s\t%s\t%s\n' "$workers" "$available" "$available" "$required" >> "$output/benchmark.tsv"
  continue
 fi
 now=$(date +%s)
 if [ "$((now-start))" -ge "$limit_seconds" ]; then
  printf '%s\ttime_budget_exhausted\t0\t0\t%s\t%s\t%s\n' "$workers" "$available" "$available" "$required" >> "$output/benchmark.tsv"
  continue
 fi
 initial=$available
 args=(--vanilla "$repo/run_deg_analysis_offline.R" --counts "$input/counts.tsv" --samples "$input/samples.tsv" --analysis "$input/analysis.tsv" --contrasts "$input/contrasts.tsv" --out "$output/workers-$workers" --workers "$workers")
 if [ -f "$input/annotation.tsv" ]; then args+=(--annotation "$input/annotation.tsv"); fi
 setsid /usr/bin/time -v -o "$output/workers-$workers.time" "$rscript" "${args[@]}" > "$output/workers-$workers.log" 2>&1 &
 active=$!
 peak=0; minimum=$available; outcome=completed; run_start=$now
 while kill -0 "$active" 2>/dev/null; do
  now=$(date +%s)
  available=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  rss=$(ps -eo pgid=,rss= | awk -v group="$active" '$1==group {sum+=$2} END {print sum+0}')
  [ "$rss" -le "$peak" ] || peak=$rss
  [ "$available" -ge "$minimum" ] || minimum=$available
  printf '%s\t%s\t%s\t%s\n' "$now" "$workers" "$rss" "$available" >> "$output/memory.tsv"
  if [ "$available" -lt "$reserve_kb" ]; then outcome=stopped_memory_guard; cleanup; break; fi
  if [ "$((now-start))" -ge "$limit_seconds" ]; then outcome=stopped_time_guard; cleanup; break; fi
  sleep 1
 done
 set +e
 wait "$active"; status=$?
 set -e
 active=
 wall=$(($(date +%s)-run_start))
 if [ "$status" -ne 0 ] && [ "$outcome" = completed ]; then outcome="failed_$status"; fi
 if [ "$workers" = 1 ]; then one_peak=$peak; fi
 if [ "$outcome" = completed ]; then
  if [ "$workers" = 1 ]; then reference=$baseline; else reference="$output/workers-1"; fi
  if ! "$rscript" --vanilla "$repo/tools/check_r_output_parity.R" "$reference" "$output/workers-$workers" > "$output/workers-$workers.parity.log" 2>&1; then outcome=parity_failed; fi
 fi
 printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$workers" "$outcome" "$wall" "$peak" "$minimum" "$initial" "$required" >> "$output/benchmark.tsv"
 # Failed baseline must not be treated as a reference for higher worker counts.
 if [ "$workers" = 1 ] && [ "$outcome" != completed ]; then break; fi
done
cat "$output/benchmark.tsv"
printf 'Artifacts: %s\n' "$output"
