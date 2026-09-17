#!/usr/bin/env bash
# Explicit native workflow entry point; input paths are relative to CONFIG.tsv.
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")" && pwd)
if [[ $# != 2 ]]; then
  echo 'Usage: run_pipeline.sh plan|local|submit CONFIG.tsv' >&2
  exit 2
fi
case "$1" in plan|local|submit) mode=$1;; *) echo 'Expected plan, local, or submit' >&2; exit 2;; esac
exec "$root/build/rnaseq" "workflow-$mode" "$2"
