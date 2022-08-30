#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Validation!" \
  -o $DATE \
  --access-level test \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/mt/986d792a448c66a8a5cfba65434e7d1ce9b1ff_1051-validation.mt \
    --header gs://cpg-validation-test/header_lines.txt
