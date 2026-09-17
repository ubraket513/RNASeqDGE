#!/usr/bin/env bash
# Stage only pinned native executables and their private shared-library closure.
set -euo pipefail
if [[ $# != 2 ]]; then echo 'Usage: stage_native_runtime.sh ALIGNMENT_PREFIX ABSENT_OUTPUT' >&2; exit 2; fi
source_prefix=$(realpath "$1")
parent=$(realpath "$(dirname "$2")")
output="$parent/$(basename "$2")"
[[ -d "$parent" && ! -e "$output" && ! -L "$output" ]] || { echo 'Output must be absent with an existing parent' >&2; exit 2; }
stage=$(mktemp -d "$parent/.native-runtime-XXXXXX")
trap 'echo "Unpublished native runtime retained: $stage" >&2' ERR
mkdir "$stage/bin" "$stage/lib" "$stage/provenance"
printf 'path\tsha256\tsource\n' > "$stage/provenance/files.tsv"
record_copy() {
  local from=$1 relative=$2 hash
  hash=$(sha256sum "$from"); hash=${hash%% *}
  if [[ -e "$stage/$relative" ]]; then
    local existing; existing=$(sha256sum "$stage/$relative"); existing=${existing%% *}
    [[ "$hash" == "$existing" ]] || { echo "Conflicting dependency $relative" >&2; return 1; }
    return
  fi
  cp -L -- "$from" "$stage/$relative"
  printf '%s\t%s\t%s\n' "$relative" "$hash" "$from" >> "$stage/provenance/files.tsv"
}
# Match the deployed P4 package pins before copying. No interpreter wrapper is
# used to query HISAT2; its small-index ELF executables are the compute interface.
[[ $("$source_prefix/bin/hisat2-align-s" --version | head -1) == *' version 2.2.1' ]]
[[ $("$source_prefix/bin/hisat2-build-s" --version | head -1) == *' version 2.2.1' ]]
[[ $("$source_prefix/bin/STAR" --version) == '2.7.10b' ]]
[[ $("$source_prefix/bin/samtools" --version | head -1) == 'samtools 1.17' ]]
[[ $("$source_prefix/bin/featureCounts" -v 2>&1 | tr -d '\r\n') == 'featureCounts v2.0.6' ]]
executables=(hisat2-align-s hisat2-build-s samtools featureCounts)
# Preserve the upstream CPU dispatcher and the ELF implementations it selects.
record_copy "$source_prefix/bin/STAR" bin/STAR
for binary in "$source_prefix"/bin/STAR-*; do
  [[ -f "$binary" && -x "$binary" ]] && executables+=("$(basename "$binary")")
done
for name in "${executables[@]}"; do
  binary="$source_prefix/bin/$name"
  [[ $(od -An -tx1 -N4 "$binary" | tr -d ' \n') == 7f454c46 ]] || { echo "Not an ELF executable: $name" >&2; exit 1; }
  record_copy "$binary" "bin/$name"
  deps=$(ldd "$binary")
  if [[ "$deps" == *'not found'* ]]; then echo "$deps" >&2; exit 1; fi
  while IFS= read -r library; do
    [[ -n "$library" ]] || continue
    resolved=$(realpath "$library")
    case "$resolved" in
      "$source_prefix"/lib/*)
        [[ $(basename "$library") != *python* ]] || { echo 'Python library in native closure' >&2; exit 1; }
        record_copy "$library" "lib/$(basename "$library")" ;;
      /lib/*|/lib64/*|/usr/lib/*|/usr/lib64/*) ;; # Host ABI libraries, recorded by ldd below.
      *) echo "Unexpected dependency outside pinned prefix or system ABI: $resolved" >&2; exit 1 ;;
    esac
  done < <(printf '%s\n' "$deps" | awk '$2 == "=>" && $3 ~ /^\// {print $3}')
  printf '%s\n' "$deps" > "$stage/provenance/$name.ldd.txt"
done
repo=$(cd "$(dirname "$0")/.." && pwd)
cp "$repo/config/alignment-linux-64.explicit.txt" "$stage/provenance/source-packages.explicit.txt"
cp "$repo/config/alignment-linux-64.provenance.json" "$stage/provenance/source-packages.provenance.json"
# The source package set includes wrapper dependencies; files.tsv is the actual
# runtime payload and contains neither Python nor the Python build wrapper.
printf 'Native compute payload; HISAT2 small indexes and uncompressed FASTQ only.\n' > "$stage/README"
(
  cd "$stage"
  find . -type f ! -name SHA256SUMS ! -name COMPLETE -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
  sha256sum -c SHA256SUMS >/dev/null
)
printf 'Python-free native runtime payload verified\n' > "$stage/COMPLETE"
mv -T --no-clobber "$stage" "$output"
[[ ! -e "$stage" ]] || { echo 'Output appeared concurrently; publication refused' >&2; exit 1; }
trap - ERR
printf 'Staged native runtime: %s\n' "$output"
