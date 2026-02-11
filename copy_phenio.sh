#!/usr/bin/env bash

curl -o phenio.db.gz https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz
gunzip phenio.db.gz
gcloud storage cp phenio.db gs://cpg-common-test/references/talos/phenio.db
