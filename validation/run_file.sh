#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Exome Validation" \
  -o $DATE \
  --access-level standard \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/exome/mt/dd7b2003026c7a6c70057a9c0f170074be6322_628-validation.mt
