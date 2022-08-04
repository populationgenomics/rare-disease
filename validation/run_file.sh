#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Validation!" \
  -o $DATE \
  --access-level test \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-validation.mt \
    --header gs://cpg-validation-test/header_lines.txt
