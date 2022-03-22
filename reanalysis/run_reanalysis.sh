#!/usr/bin/env bash

set -ex

# set the date, or provide a default
OUTPUT_DATE_FOLDER=${1:-"2021-09-03"}

# run
analysis-runner \
  --dataset acute-care \
  --description "run reanalysis rewrite" \
  -o "reanalysis/${OUTPUT_DATE_FOLDER}" \
  --access-level test \
  reanalysis/reanalysis_wrapper.py \
    --conf gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json \
    --matrix gs://cpg-acute-care-main/mt/acute-care.mt \
    --panel_genes gs://cpg-acute-care-test/reanalysis/pre_panelapp_mendeliome.json \
    --ped gs://cpg-acute-care-test/reanalysis/acute_care_singleton_pedigree.json
