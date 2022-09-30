#!/usr/bin/env bash

set -ex

DATE=${1:-$(date +%F)}

analysis-runner \
  --dataset validation \
  --description "Run Exome Validation!" \
  -o $DATE \
  --access-level test \
  validation/validation_runner.py \
    -i gs://cpg-validation-test/validation/copy/copy_of_exome_628.mt \
    -s gs://cpg-validation-test/GRCh38_regions/twist/definition.tsv
