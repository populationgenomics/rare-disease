#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Genome Validation" \
  -o $DATE \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows \
  --access-level standard \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/mt/6ae7ac744240e459f9e38f794631957c066ae4_1359-validation.mt \
    -s gs://cpg-validation-test/GRCh38_regions
#    --dry_run
