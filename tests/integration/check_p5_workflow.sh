#!/usr/bin/env bash
set -euo pipefail
exe=$(realpath "${1:-build/rnaseq}")
repo=$(pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/bin" "$tmp/scheduler"
cp -R tests/fixtures/p0 "$tmp/p0"
cp tests/fixtures/p2/genes.tsv "$tmp/genes.tsv"
cp run_deg_analysis_offline.R "$tmp/run_deg_analysis_offline.R"
cat > "$tmp/bin/tool" <<'STUB'
#!/bin/sh
name=${0##*/}
if [ "$1" = --version ] || [ "$1" = -v ]; then
 case "$name" in hisat2*) echo "$name version 2.2.1";; samtools) echo 'samtools 1.17';; featureCounts) echo 'featureCounts v2.0.6';; esac
 exit 0
fi
printf '%s\n' "$name" >> "$EVENTS"
if [ "$name" = hisat2-build-s ]; then
 if [ -n "${BLOCK:-}" ]; then
  trap '' INT TERM
  sleep 200 &
  printf '%s %s\n' "$$" "$!" > "$BLOCK"
  wait
 fi
 for prefix; do :; done
 for i in 1 2 3 4 5 6 7 8; do echo index > "$prefix.$i.ht2"; done
elif [ "$name" = hisat2-align-s ]; then
 if [ -n "${PARALLEL_DIR:-}" ]; then
  printf 'start %s\n' "$$" >> "$PARALLEL_DIR/events"
  if [ "${PARALLEL_MODE:-overlap}" = overlap ]; then
   sleep .2
  else
   trap '' INT TERM
   setsid sleep 200 &
   printf '%s %s\n' "$$" "$!" > "$PARALLEL_DIR/worker.$$"
   if [ "$PARALLEL_MODE" = failure ]; then
    while [ "$(find "$PARALLEL_DIR" -name 'worker.*' | wc -l)" -lt 2 ]; do sleep .01; done
    if mkdir "$PARALLEL_DIR/first-failure" 2>/dev/null; then exit 19; fi
   fi
   wait
  fi
  printf 'end %s\n' "$$" >> "$PARALLEL_DIR/events"
 fi
 while [ "$#" -gt 0 ]; do if [ "$1" = -S ]; then shift; echo sam > "$1"; fi; shift; done
elif [ "$name" = samtools ]; then
 [ "$1" = quickcheck ] && exit 0
 while [ "$#" -gt 0 ]; do if [ "$1" = -o ]; then shift; echo bam > "$1"; fi; shift; done
elif [ "$name" = featureCounts ]; then
 while [ "$#" -gt 0 ]; do if [ "$1" = -o ]; then shift; out=$1; fi; bam=$1; shift; done
 printf '# stub fixture\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t%s\n' "$bam" > "$out"
 for gene in gene_A.1 gene_B.2 gene_C gene_zero; do printf '%s\tchr1\t1\t10\t+\t10\t1\n' "$gene" >> "$out"; done
 printf 'Status\t%s\nAssigned\t4\nUnassigned_Unmapped\t0\n' "$bam" > "$out.summary"
fi
STUB
chmod +x "$tmp/bin/tool"
for tool in hisat2-build-s hisat2-align-s samtools featureCounts; do cp "$tmp/bin/tool" "$tmp/bin/$tool"; done
cat > "$tmp/Rscript" <<'STUB'
#!/bin/sh
[ "${OMP_NUM_THREADS:-}" = 1 ] && [ "${OPENBLAS_NUM_THREADS:-}" = 1 ] || exit 71
[ -z "${R_HOME:-}" ] && [ -z "${R_LIBS:-}" ] && [ -z "${LD_LIBRARY_PATH:-}" ] || exit 72
[ "$R_LIBS_USER" = /dev/null ] && [ "$R_LIBS_SITE" = /dev/null ] || exit 73
case "$2" in */r-preflight.R) exit 0;; esac
printf 'R\n' >> "$EVENTS"
[ -z "${R_FAIL:-}" ] || exit 19
while [ "$#" -gt 0 ]; do if [ "$1" = --out ]; then shift; out=$1; fi; shift; done
mkdir "$out"
printf 'stub boundary only\n' > "$out/results.tsv"
STUB
chmod +x "$tmp/Rscript"
printf 'pinned fixture lock\n' > "$tmp/lock"
export EVENTS="$tmp/events"
config() {
 local output=$1 run=$2
 printf 'key\tvalue\nsamples\t%s/p0/samples.tsv\nruns\t%s/p0/runs.tsv\nreferences\t%s/p0/references.tsv\nanalysis\t%s/p0/analysis.tsv\ncontrasts\t%s/p0/contrasts.tsv\ngenes\t%s/genes.tsv\nbin_dir\t%s/bin\nrscript\t%s/Rscript\nr_script\t%s/run_deg_analysis_offline.R\ntool_lock\t%s/lock\nr_lock\t%s/lock\nrun_dir\t%s\n' "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$tmp" "$run" > "$output"
}
reject() { if "$@" > "$tmp/reject.log" 2>&1; then echo "unexpected success: $*" >&2; exit 1; fi; }
# Submission describes the target allocation, not the login host's affinity.
config "$tmp/affinity.tsv" "$tmp/affinity-run"
mkdir "$tmp/affinity-scheduler"
printf '#!/bin/sh\nprintf "101\\n"\n' > "$tmp/affinity-scheduler/sbatch"
printf '#!/bin/sh\nexit 0\n' > "$tmp/affinity-scheduler/scancel"
chmod +x "$tmp/affinity-scheduler/"*
printf 'slurm_bin_dir\t%s/affinity-scheduler\nslurm_cpus\t32\nslurm_mem_mb\t4096\nslurm_time\t00:30:00\nthreads\t8\nlocal_cpus\t1\n' "$tmp" >> "$tmp/affinity.tsv"
login_cpu=$(awk '/Cpus_allowed_list/ {split($2,a,/[,-]/); print a[1]}' /proc/self/status)
taskset -c "$login_cpu" "$exe" workflow-plan "$tmp/affinity.tsv" > "$tmp/affinity.plan"
grep -q 'resource_scope.*slurm_requested' "$tmp/affinity.plan"
taskset -c "$login_cpu" "$exe" workflow-submit "$tmp/affinity.tsv" > "$tmp/affinity.log"
grep -qx -- '--cpus-per-task=32' "$tmp/affinity-run/submission/index.argv"
reject taskset -c "$login_cpu" "$exe" workflow-local "$tmp/affinity.tsv"
grep -q 'CPU budget' "$tmp/reject.log"
reject env SLURM_CPUS_PER_TASK=32 taskset -c "$login_cpu" "$exe" workflow-task "$tmp/affinity.tsv" index
grep -q 'CPU budget' "$tmp/reject.log"
sed 's/slurm_cpus\t32/slurm_cpus\t4/' "$tmp/affinity.tsv" > "$tmp/affinity-small.tsv"
reject taskset -c "$login_cpu" "$exe" workflow-plan "$tmp/affinity-small.tsv"
grep -q 'slurm_cpus smaller than threads/workers' "$tmp/reject.log"
reject taskset -c "$login_cpu" "$exe" workflow-submit "$tmp/affinity-small.tsv"
grep -q 'slurm_cpus smaller than threads/workers' "$tmp/reject.log"
printf 'PASS P5 submission affinity: target Slurm allocation independent of login CPU limits\n'
if test "${P5_AFFINITY_ONLY:-0}" = 1; then exit 0; fi
config "$tmp/config.tsv" "$tmp/run ; literal"
reject "$exe" workflow-task "$tmp/config.tsv"
reject "$exe" workflow-unknown "$tmp/config.tsv"
"$exe" workflow-plan "$tmp/config.tsv" > "$tmp/plan"
test ! -e "$tmp/events" && test ! -e "$tmp/run ; literal"
grep -q 'threads.*1' "$tmp/plan"
# Local scheduling resolves CPU and reservation budgets before execution.
config "$tmp/budget.tsv" "$tmp/budget-run"
printf 'local_jobs\t4\nlocal_cpus\t4\nlocal_mem_mb\t200\nlocal_job_mem_mb\t100\n' >> "$tmp/budget.tsv"
SLURM_CPUS_PER_TASK=4 "$exe" workflow-plan "$tmp/budget.tsv" > "$tmp/budget.plan"
grep -q 'local_resolved_jobs.*2' "$tmp/budget.plan"
SLURM_CPUS_PER_TASK=1 "$exe" workflow-plan "$tmp/budget.tsv" > "$tmp/budget-one.plan"
grep -q 'local_resolved_jobs.*1' "$tmp/budget-one.plan"
SLURM_CPUS_PER_TASK=2 SLURM_MEM_PER_CPU=50 "$exe" workflow-plan "$tmp/budget.tsv" > "$tmp/budget-memory.plan"
grep -q 'local_resolved_jobs.*1' "$tmp/budget-memory.plan"
config "$tmp/no-mem.tsv" "$tmp/no-mem"
printf 'local_jobs\t4\nlocal_cpus\t4\n' >> "$tmp/no-mem.tsv"
"$exe" workflow-plan "$tmp/no-mem.tsv" > "$tmp/no-mem.plan"
grep -q 'local_resolved_jobs.*1' "$tmp/no-mem.plan"
grep -q 'serial fallback' "$tmp/no-mem.plan"
config "$tmp/too-small.tsv" "$tmp/too-small"
printf 'local_mem_mb\t50\nlocal_job_mem_mb\t100\n' >> "$tmp/too-small.tsv"
reject "$exe" workflow-plan "$tmp/too-small.tsv"
grep -q 'local_job_mem_mb exceeds' "$tmp/reject.log"
reject env SLURM_CPUS_PER_TASK=0 "$exe" workflow-plan "$tmp/config.tsv"
R_LIBS=/unwanted R_LIBS_USER=/unwanted R_LIBS_SITE=/unwanted R_HOME=/unwanted LD_LIBRARY_PATH=/unwanted "$exe" workflow-local "$tmp/config.tsv" > "$tmp/local.log"
run="$tmp/run ; literal"
test -f "$run/analysis/STAGE.tsv"
test -d "$run/profiles"
grep -q executed "$run/profiles/"*.tsv
grep -q '.staging-.*aligned.bam' "$run/merge/inputs.tsv"
cp "$EVENTS" "$tmp/prior-events"
"$exe" workflow-local "$tmp/config.tsv" > "$tmp/resume.log"
cmp "$EVENTS" "$tmp/prior-events"
grep -q reused "$run/profiles/"*.tsv
printf x >> "$tmp/p0/reads/single.fastq"
reject "$exe" workflow-local "$tmp/config.tsv"
grep -q 'new run directory' "$tmp/reject.log"
cp "$repo/tests/fixtures/p0/reads/single.fastq" "$tmp/p0/reads/single.fastq"
printf '\n' >> "$tmp/bin/hisat2-align-s"
reject "$exe" workflow-local "$tmp/config.tsv"
cp "$tmp/bin/tool" "$tmp/bin/hisat2-align-s"
printf 'threads\t2\n' >> "$tmp/config.tsv"
reject "$exe" workflow-local "$tmp/config.tsv"
sed -i '$d' "$tmp/config.tsv"
printf corrupt >> "$run/align/run_t2/counts.txt"
"$exe" workflow-local "$tmp/config.tsv" > "$tmp/recover.log"
grep -q quarantine "$tmp/recover.log"
test "$(grep -c '^R$' "$EVENTS")" = 2
rm "$run/merge/STAGE.tsv"
"$exe" workflow-local "$tmp/config.tsv" > "$tmp/partial.log"
grep -q 'merge.invalid' "$tmp/partial.log"
test "$(grep -c '^R$' "$EVENTS")" = 3
# Two alignment workers overlap, obey reservations, merge deterministically and resume.
mkdir "$tmp/parallel-tools"
PARALLEL_DIR="$tmp/parallel-tools" "$exe" workflow-local "$tmp/budget.tsv" > "$tmp/parallel.log"
awk '$1=="start" {n++; if(n>max)max=n} $1=="end" {n--} END {exit !(max==2 && n==0)}' "$tmp/parallel-tools/events"
cmp "$tmp/budget-run/merge/counts.tsv" "$run/merge/counts.tsv"
cp "$tmp/parallel-tools/events" "$tmp/parallel-events"
PARALLEL_DIR="$tmp/parallel-tools" "$exe" workflow-local "$tmp/budget.tsv" > "$tmp/parallel-resume.log"
cmp "$tmp/parallel-tools/events" "$tmp/parallel-events"
# Failure and parent signals kill/reap descendants even after setsid().
for mode in failure signal signal-int; do
 mkdir "$tmp/$mode-tools"
 config "$tmp/$mode-parallel.tsv" "$tmp/$mode-parallel"
 printf 'local_jobs\t2\nlocal_cpus\t2\nlocal_mem_mb\t200\nlocal_job_mem_mb\t100\n' >> "$tmp/$mode-parallel.tsv"
 PARALLEL_DIR="$tmp/$mode-tools" PARALLEL_MODE="$mode" "$exe" workflow-local "$tmp/$mode-parallel.tsv" > "$tmp/$mode-parallel.log" 2>&1 &
 concurrent=$!
 for _ in $(seq 1 500); do
  test "$(find "$tmp/$mode-tools" -name 'worker.*' | wc -l)" -ge 2 && break
  sleep .01
 done
 test "$(find "$tmp/$mode-tools" -name 'worker.*' | wc -l)" -ge 2
 if test "$mode" = signal; then kill -TERM "$concurrent"; fi
 if test "$mode" = signal-int; then kill -INT "$concurrent"; fi
 set +e
 wait "$concurrent"; status=$?
 set -e
 if test "$mode" = signal; then test "$status" = 143; elif test "$mode" = signal-int; then test "$status" = 130; else test "$status" = 19; fi
 for record in "$tmp/$mode-tools"/worker.*; do
  read -r tool descendant < "$record"
  test ! -e "/proc/$tool"
  test ! -e "/proc/$descendant"
 done
 test ! -e "$tmp/$mode-parallel/merge/STAGE.tsv"
 test ! -e "$tmp/$mode-parallel/analysis/STAGE.tsv"
done
# Failure diagnostics survive; restart retries only the failed R boundary.
config "$tmp/fail.tsv" "$tmp/fail"
reject env R_FAIL=1 "$exe" workflow-local "$tmp/fail.tsv"
test ! -f "$tmp/fail/analysis/STAGE.tsv"
"$exe" workflow-local "$tmp/fail.tsv" > "$tmp/restart.log"
test -f "$tmp/fail/analysis/STAGE.tsv"
# Shared cache reuse and simultaneous run/cache protection.
config "$tmp/cache1.tsv" "$tmp/cache1-run"
printf 'index_cache\t%s/cache\n' "$tmp" >> "$tmp/cache1.tsv"
"$exe" workflow-local "$tmp/cache1.tsv" > "$tmp/cache1.log"
config "$tmp/cache2.tsv" "$tmp/cache2-run"
printf 'index_cache\t%s/cache\n' "$tmp" >> "$tmp/cache2.tsv"
"$exe" workflow-local "$tmp/cache2.tsv" > "$tmp/cache2.log"
grep -q 'reuse.*cache/' "$tmp/cache2.log"
# Descendant cancellation is exercised through native workflow -> P4 -> tool.
config "$tmp/block.tsv" "$tmp/block-run"
printf 'index_cache\t%s/block-cache\n' "$tmp" >> "$tmp/block.tsv"
BLOCK="$tmp/blocked" "$exe" workflow-local "$tmp/block.tsv" > "$tmp/block.log" 2>&1 &
runner=$!
for _ in $(seq 1 300); do test -f "$tmp/blocked" && break; sleep .01; done
test -f "$tmp/blocked"
reject "$exe" workflow-local "$tmp/block.tsv"
grep -q 'busy lock' "$tmp/reject.log"
config "$tmp/contend.tsv" "$tmp/contend-run"
printf 'index_cache\t%s/block-cache\n' "$tmp" >> "$tmp/contend.tsv"
reject "$exe" workflow-local "$tmp/contend.tsv"
grep -q 'busy lock.*block-cache' "$tmp/reject.log"
read -r child grandchild < "$tmp/blocked"
kill -TERM "$runner"
set +e
wait "$runner"; status=$?
set -e
test "$status" = 143
! kill -0 "$child" 2>/dev/null
# Grandchild can briefly remain a zombie adopted by PID 1, but must not run.
if test -e "/proc/$grandchild/stat"; then test "$(awk '{print $3}' "/proc/$grandchild/stat")" = Z; fi
test "$(find "$tmp/block-cache" -name STAGE.tsv | wc -l)" = 0
# Fake scheduler records exact argv and simulates failures after an accepted job.
cat > "$tmp/scheduler/sbatch" <<'STUB'
#!/bin/sh
n=0; [ ! -f "$SCHED/count" ] || n=$(cat "$SCHED/count")
n=$((n+1)); echo "$n" > "$SCHED/count"
printf '%s\n' "$@" > "$SCHED/args-$n"
[ "${SUBMIT_FAIL:-0}" != "$n" ] || exit 23
echo 'scheduler informational warning' >&2
if [ "${EARLY_TASK:-0}" = 1 ] && [ "$n" = 1 ]; then
 for script; do :; done
 /bin/sh "$script" > "$SCHED/early-task.log" 2>&1 &
 echo "$!" > "$SCHED/early-task.pid"
fi
echo "$((100+n))"
STUB
cat > "$tmp/scheduler/scancel" <<'STUB'
#!/bin/sh
printf '%s\n' "$@" >> "$SCHED/cancelled"
STUB
chmod +x "$tmp/scheduler/"*
export SCHED="$tmp/scheduler"
config "$tmp/submit.tsv" "$tmp/submit"
printf 'slurm_bin_dir\t%s/scheduler\nslurm_cpus\t2\nslurm_mem_mb\t4096\nslurm_time\t00:30:00\nslurm_concurrency\t2\nslurm_partition\tdebug\n' "$tmp" >> "$tmp/submit.tsv"
EARLY_TASK=1 "$exe" workflow-submit "$tmp/submit.tsv" > "$tmp/submit.log"
for _ in $(seq 1 500); do test -f "$tmp/submit/index/STAGE.tsv" && break; sleep .01; done
test -f "$tmp/submit/index/STAGE.tsv"
grep -q 'scheduler informational warning' "$tmp/submit/submission/index.submit.stderr"
! grep -q 'busy lock' "$SCHED/early-task.log"
grep -qx -- '--parsable' "$SCHED/args-1"
grep -qx -- '--dependency=afterok:101' "$SCHED/args-2"
grep -qx -- '--array=0-4%2' "$SCHED/args-2"
grep -qx -- '--dependency=afterok:102' "$SCHED/args-3"
grep -qx -- '--cpus-per-task=2' "$SCHED/args-1"
grep -qx -- '--mem=4096' "$SCHED/args-1"
grep -qx -- '--time=00:30:00' "$SCHED/args-1"
reject "$exe" workflow-submit "$tmp/submit.tsv"
test "$(cat "$SCHED/count")" = 3
# Execute generated scripts locally with scheduler variables to test mapping.
bash "$tmp/submit/submission/index.sh" > "$tmp/task-index.log"
SLURM_ARRAY_TASK_ID=0 bash "$tmp/submit/submission/array.sh" > "$tmp/task-0.log" &
task0=$!
SLURM_ARRAY_TASK_ID=1 bash "$tmp/submit/submission/array.sh" > "$tmp/task-1.log" &
task1=$!
wait "$task0"
wait "$task1"
for i in 2 3 4; do SLURM_ARRAY_TASK_ID=$i bash "$tmp/submit/submission/array.sh" > "$tmp/task-$i.log"; done
bash "$tmp/submit/submission/finish.sh" > "$tmp/task-finish.log"
test -f "$tmp/submit/analysis/STAGE.tsv"
sed "s|$tmp/submit$|$tmp/rollback|" "$tmp/submit.tsv" > "$tmp/rollback.tsv"
reject env SUBMIT_FAIL=5 "$exe" workflow-submit "$tmp/rollback.tsv"
grep -qx 104 "$SCHED/cancelled"
test -f "$tmp/rollback/submission/FAILED"
printf 'PASS P5 workflow: DAG/resume, source/config/tool invalidation, corrupt/partial recovery, restart, cache reuse, locks, descendant cancellation, Slurm dependencies/resources/tasks/rollback\n'

# A cancelled preflight hash must not leave a large-file checksum reader alive.
mkdir "$tmp/hash-bin"
cat > "$tmp/hash-bin/sha256sum" <<'STUB'
#!/bin/sh
printf '%s\n' "$$" > "$HASH_PID"
while :; do :; done
STUB
chmod +x "$tmp/hash-bin/sha256sum"
HASH_PID="$tmp/hash.pid" PATH="$tmp/hash-bin:$PATH" "$exe" workflow-plan "$tmp/config.tsv" > "$tmp/hash.log" 2>&1 &
hash_driver=$!
for _ in $(seq 1 200); do test -f "$tmp/hash.pid" && break; sleep .01; done
test -f "$tmp/hash.pid"
hash_child=$(cat "$tmp/hash.pid")
kill -TERM "$hash_driver"
wait "$hash_driver" 2>/dev/null || true
for _ in $(seq 1 100); do
 test ! -e "/proc/$hash_child/stat" && break
 test "$(awk '{print $3}' "/proc/$hash_child/stat")" = Z && break
 sleep .01
done
if test -e "/proc/$hash_child/stat"; then test "$(awk '{print $3}' "/proc/$hash_child/stat")" = Z; fi
printf 'PASS P5 hash cancellation: checksum child stopped with driver\n'
