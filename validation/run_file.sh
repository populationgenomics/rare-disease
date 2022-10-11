#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Genome Validation" \
  -o $DATE \
  --access-level standard \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/mt/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-validation.mt \
    -s gs://cpg-validation-test/GRCh38_regions \
    --no_post
