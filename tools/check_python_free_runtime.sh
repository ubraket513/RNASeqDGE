#!/usr/bin/env bash
# Check the deployed prefixes, not historical source-package provenance.
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: check_python_free_runtime.sh R_PREFIX NATIVE_PREFIX' >&2; exit 2; }
r_prefix=$(realpath "$1")
native_prefix=$(realpath "$2")
for prefix in "$r_prefix" "$native_prefix"; do
  [[ -d "$prefix/bin" ]] || exit 1
  bad=$(find "$prefix" \( -type f -o -type l \) \( -name 'libpython*' -o \( -executable \( -name 'python*' -o -name 'pypy*' \) \) \) -print)
  [[ -z "$bad" ]] || { printf 'Python payload found:\n%s\n' "$bad" >&2; exit 1; }
done
[[ ! -e "$native_prefix/bin/hisat2-build" && ! -e "$native_prefix/bin/hisat2" ]]
(cd "$native_prefix" && sha256sum --quiet -c SHA256SUMS)
for binary in "$native_prefix"/bin/*; do
  if [[ $(od -An -tx1 -N4 "$binary" | tr -d ' \n') == 7f454c46 ]]; then
    libraries=$(ldd "$binary")
    [[ "$libraries" != *'not found'* && "$libraries" != *libpython* ]] || exit 1
    while read -r library; do
      case "$library" in "$native_prefix"/lib/*|"$native_prefix"/bin/../lib/*|/lib/*|/lib64/*|/usr/lib/*|/usr/lib64/*) ;; *) echo "Unexpected runtime dependency: $library" >&2; exit 1;; esac
    done < <(printf '%s\n' "$libraries" | awk '$2 == "=>" && $3 ~ /^\// {print $3}')
  fi
done
printf 'PASS: R and native prefixes contain no Python interpreter/library; native closure and hashes verified\n'
