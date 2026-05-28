#!/usr/bin/env bash
# check_gvcfs_from_files.sh <file_paths.txt>

set -uo pipefail
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

check_one() {
  local tbi_file="$1"

  # Skip empty lines
  [[ -z "$tbi_file" ]] && return

  local gvcf="${tbi_file%.tbi}"
  local tmpdir stats last creation_date

  tmpdir=$(mktemp -d)
  trap 'rm -rf "$tmpdir"' RETURN

  # Fetch the creation date using gsutil stat
  creation_date=$(gsutil stat "${tbi_file}" 2>/dev/null | grep "Creation time:" | sed 's/^[[:space:]]*Creation time:[[:space:]]*//')
  # Fallback if date can't be fetched
  [[ -z "$creation_date" ]] && creation_date="UNKNOWN_DATE"

  stats=$(cd "$tmpdir" && bcftools index --stats "${tbi_file}" 2>/dev/null) || {
    echo "ERROR ${gvcf} last=N/A date=${creation_date}"
    return
  }

  # Extract the last contig (config) for ALL files
  last=$(tail -n1 <<<"$stats" | cut -f1)
  [[ -z "$last" ]] && last="NONE"

  if ! grep -q $'^chrX\t' <<<"$stats"; then
    echo "TRUNCATED ${gvcf}  last=${last}  date=${creation_date}"
  else
    echo "OK  ${gvcf}  last=${last}  date=${creation_date}"
  fi
}
export -f check_one

# Check if file argument is provided
if [[ -z "${1:-}" || ! -f "$1" ]]; then
  echo "Usage: $0 <file_paths.txt>"
  exit 1
fi

# Read individual file paths directly, ignore empty lines, and process in parallel
grep -v '^$' "$1" | xargs -P 16 -I{} bash -c 'check_one "$@"' _ {}