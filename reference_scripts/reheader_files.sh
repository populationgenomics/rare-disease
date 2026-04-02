#!/usr/bin/env bash

IN=$1
OUT=$2

gcloud storage cat "${IN}" | gunzip | awk '{if(NR==1){print "#CHROM\tPOS\tREF\tALT\tavis\tphred"}else{print $0}}' | bgzip -c | gcloud storage cp "${OUT}"
