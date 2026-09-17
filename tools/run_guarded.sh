#!/usr/bin/env bash
# Run one argv command in a fresh session and forward cancellation to its group.
set -uo pipefail
[[ $# -gt 0 ]] || exit 2
setsid "$@" &
child=$!
cancel() {
  local signal=$1 code=$2 watchdog
  trap '' INT TERM
  kill -"$signal" -- -"$child" 2>/dev/null || true
  (sleep 15; kill -KILL -- -"$child" 2>/dev/null || true) &
  watchdog=$!
  wait "$child" 2>/dev/null || true
  kill "$watchdog" 2>/dev/null || true
  wait "$watchdog" 2>/dev/null || true
  exit "$code"
}
trap 'cancel INT 130' INT
trap 'cancel TERM 143' TERM
wait "$child"
exit $?
