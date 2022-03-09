#!/usr/bin/env bash

# takes a file which was pre-processed by Hail to include a compound CSQ field
# use this script to add a required Description to the file
file_in=$1
file_ext="${file_in#*.}"
filename="${file_in%%.*}"
file_out="${filename}_reheadered.${file_ext}"

bcftools view -h "${file_in}" | \
  sed 's/##INFO=<ID=COMPOUND_CSQ,Number=1,Type=String,Description="">/##INFO=<ID=COMPOUND_CSQ,Number=1,Type=String,Description="Format: 'Gene\|Transcript\|Consequence'">/' > new_header

bcftools reheader -h new_header --threads 4 -o "${file_out}" "${file_in}"

tabix "${file_out}"
