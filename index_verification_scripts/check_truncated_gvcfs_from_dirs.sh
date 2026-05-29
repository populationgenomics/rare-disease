#!/usr/bin/env bash
# check_truncated_gvcfs_from_dirs.sh <paths_file>

set -uo pipefail
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

check_one() {
  local tbi_file="$1"
  local gvcf="${tbi_file%.tbi}"
  local tmpdir stats

  tmpdir=$(mktemp -d)
  trap 'rm -rf "$tmpdir"' RETURN

  stats=$(cd "$tmpdir" && bcftools index --stats "${tbi_file}" 2>/dev/null) || {
    echo "ERROR	${gvcf}"; return
  }

  if ! grep -q $'^chrX\t' <<<"$stats"; then
    last=$(tail -n1 <<<"$stats" | cut -f1)
    echo "TRUNCATED	${gvcf}	last=${last}"
  else
    echo "OK	${gvcf}"
  fi
}
export -f check_one

# Read directories from the input file, fetch .tbi files, and process them in parallel
while IFS= read -r bucket_path || [[ -n "$bucket_path" ]]; do
  if [[ -n "$bucket_path" ]]; then
    gsutil ls "${bucket_path%%/}/*.tbi" 2>/dev/null
  fi
done < "$1" | xargs -P 16 -I{} bash -c 'check_one "$@"' _ {}
