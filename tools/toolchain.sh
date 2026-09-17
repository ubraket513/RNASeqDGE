#!/usr/bin/env bash
# Bootstrap explicit environments without requiring Python or an existing R.
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
if [[ ${1:-} != stage ]]; then exec Rscript --vanilla "$root/tools/toolchain.R" "$@"; fi
shift
kind=${1:?Usage: toolchain.sh stage runtime-r|alignment|p0 --prefix PATH [--mamba EXECUTABLE]}
shift
case "$kind" in runtime-r|alignment|p0) ;; *) echo 'unknown environment kind' >&2; exit 2;; esac
prefix= manager=mamba
while (($#)); do
  case "$1" in
    --prefix) [[ -z "$prefix" && $# -ge 2 ]] || exit 2; prefix=$2; shift 2;;
    --mamba) [[ $# -ge 2 ]] || exit 2; manager=$2; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[[ -n "$prefix" && ! -e "$prefix" && ! -L "$prefix" && -d "$(dirname -- "$prefix")" ]] || { echo 'prefix must be absent with existing parent' >&2; exit 2; }
lock="$root/config/$kind-linux-64.explicit.txt"
grep -qx '@EXPLICIT' "$lock"
# Explicit locks select package archives; this command never invokes a solver.
env -u LD_LIBRARY_PATH -u R_HOME -u R_LIBS R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$manager" create --yes --prefix "$prefix" --file "$lock"
rscript=Rscript
[[ ! -x "$prefix/bin/Rscript" ]] || rscript="$prefix/bin/Rscript"
exec "$rscript" --vanilla "$root/tools/toolchain.R" verify --kind "$kind" --prefix "$prefix"
