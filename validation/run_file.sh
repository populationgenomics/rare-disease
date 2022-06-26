#!/usr/bin/env bash

set -ex

analysis-runner \
  --dataset validation \
  --description "Convert VCF" \
  -o "2022-06-24" \
  --access-level test \
  validation/validation_runner.py \
    -i gs://cpg-validation-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-validation.mt \
    -h gs://cpg-validation-test/header_lines.txt
