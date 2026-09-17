#!/usr/bin/env bash
# Runs an explicitly staged reduced-data workflow; never downloads or selects samples.
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
[[ $# == 1 ]] || { echo 'Usage: tests/test_subsampled.sh CONFIG.tsv' >&2; exit 2; }
exec bash "$root/run_pipeline.sh" local "$1"
