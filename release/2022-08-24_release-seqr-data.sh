#!/usr/bin/env bash
# https://github.com/populationgenomics/seqr-private/issues/8

# 121.64 GiB (~ $36 egress)
gsutil -m -u perth-neuro cp \
	'gs://cpg-perth-neuro-main/cram/CPG13136.cram' 'gs://cpg-perth-neuro-main/cram/CPG13136.cram.crai' \
	'gs://cpg-perth-neuro-main/cram/CPG13334.cram' 'gs://cpg-perth-neuro-main/cram/CPG13334.cram.crai' \
	'gs://cpg-perth-neuro-main/cram/CPG13342.cram' 'gs://cpg-perth-neuro-main/cram/CPG13342.cram.crai' \
	'gs://cpg-perth-neuro-main/cram/CPG222430.cram' 'gs://cpg-perth-neuro-main/cram/CPG222430.cram.crai' \
	'gs://cpg-perth-neuro-main/cram/CPG222448.cram' 'gs://cpg-perth-neuro-main/cram/CPG222448.cram.crai' \
	'gs://cpg-perth-neuro-main/cram/CPG222455.cram' 'gs://cpg-perth-neuro-main/cram/CPG222455.cram.crai' \
	'gs://cpg-perth-neuro-release/'

# 16.7 GiB (~ $5 AUD egress)
gsutil -m -u ravenscroft-rdstudy cp \
	'gs://cpg-ravenscroft-rdstudy-main/cram/CPG202655.cram' 'gs://cpg-ravenscroft-rdstudy-main/cram/CPG202655.cram.crai' \
	'gs://cpg-ravenscroft-rdstudy-release/'

# 816.44 GiB (~ $240 AUD egress)
gsutil -m -u mito-disease cp \
	'gs://cpg-mito-disease-main/cram/CPG210815.cram' 'gs://cpg-mito-disease-main/cram/CPG210815.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211409.cram' 'gs://cpg-mito-disease-main/cram/CPG211409.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211359.cram' 'gs://cpg-mito-disease-main/cram/CPG211359.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211284.cram' 'gs://cpg-mito-disease-main/cram/CPG211284.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210849.cram' 'gs://cpg-mito-disease-main/cram/CPG210849.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211375.cram' 'gs://cpg-mito-disease-main/cram/CPG211375.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211334.cram' 'gs://cpg-mito-disease-main/cram/CPG211334.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210807.cram' 'gs://cpg-mito-disease-main/cram/CPG210807.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211227.cram' 'gs://cpg-mito-disease-main/cram/CPG211227.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211326.cram' 'gs://cpg-mito-disease-main/cram/CPG211326.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211276.cram' 'gs://cpg-mito-disease-main/cram/CPG211276.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210864.cram' 'gs://cpg-mito-disease-main/cram/CPG210864.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210872.cram' 'gs://cpg-mito-disease-main/cram/CPG210872.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211268.cram' 'gs://cpg-mito-disease-main/cram/CPG211268.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211383.cram' 'gs://cpg-mito-disease-main/cram/CPG211383.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211243.cram' 'gs://cpg-mito-disease-main/cram/CPG211243.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210856.cram' 'gs://cpg-mito-disease-main/cram/CPG210856.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211300.cram' 'gs://cpg-mito-disease-main/cram/CPG211300.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211391.cram' 'gs://cpg-mito-disease-main/cram/CPG211391.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211292.cram' 'gs://cpg-mito-disease-main/cram/CPG211292.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211367.cram' 'gs://cpg-mito-disease-main/cram/CPG211367.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211250.cram' 'gs://cpg-mito-disease-main/cram/CPG211250.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210799.cram' 'gs://cpg-mito-disease-main/cram/CPG210799.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211342.cram' 'gs://cpg-mito-disease-main/cram/CPG211342.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210831.cram' 'gs://cpg-mito-disease-main/cram/CPG210831.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211318.cram' 'gs://cpg-mito-disease-main/cram/CPG211318.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG210823.cram' 'gs://cpg-mito-disease-main/cram/CPG210823.cram.crai' \
	'gs://cpg-mito-disease-main/cram/CPG211235.cram' 'gs://cpg-mito-disease-main/cram/CPG211235.cram.crai' \
	'gs://cpg-mito-disease-release/'
