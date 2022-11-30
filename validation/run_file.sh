#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Genome Validation" \
  -o $DATE \
  --config validation/validation_conf.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
  --access-level test \
  validation/validation_runner.py \
    -i gs://cpg-validation-test/validation/copy/copy_of_1359.mt \
    -s gs://cpg-validation-test/GRCh38_regions \
    --dry_run
