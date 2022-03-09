# Design Doc

## Aim

Take in variant data, and run a variant prioritisation algorithm. Highlight a minimal set of variants per sample, 
based on strict criteria

In general this is focused on True Positives, and is happy to allow False Negatives, making a highly specific 
actionable dataset, at the potential expense of sensitivity

This is currently in an MVP phase, where the core logic should be accurate but the input and output are only 
rough drafts

---

## Concept

This product runs as a single Hail Batch workflow, with a number of Sub-Batches being generated 
for individual containerised steps. The intended runtime environment is the 
[analysis-runner](https://github.com/populationgenomics/analysis-runner), which will handle authentication.

1. Query PanelApp to obtain gene lists & MOI data
2. Load variant data in Hail Query, annotate, filter, and assign categories of interest
3. Use Slivar to identify and annotate compound-heterozygous variants within the joint-call
4. Combine categorised data, gene data, and compound hets to variants of interest

This requires pedigree data to describe the family structures. Initially this will be provided using a file, but 
code is included that can query the [sample metadata](https://github.com/populationgenomics/sample-metadata) API.

Input required:
1. Variant Data in either Hail MatrixTable or VCF format
2. Pedigree
3. Configuration file (default provided in this repo)


Outputs generated:
1. PanelApp data
2. MatrixTable containing annotated, unfiltered data
3. VCF containing annotated categorised variants
4. VCF containing compound-heterozygous variants
5. HTML report of confirmed categorised variants

The outputs will be written to the 
[CPG GCP bucket](https://github.com/populationgenomics/team-docs/tree/main/storage_policies) corresponding 
to the cohort being analysed. The intermediate files and outputs will be stored in the project's core bucket, 
with the report HTML additionally available in the `-web` suffixed bucket 

---

# Instructions

Run using a call to the analysis-runner, e.g.

```bash
PanelApp_Date=${1:-"2021-09-03"}
analysis-runner \
  --dataset acute-care \
  --description "run reanalysis draft" \
  -o "reanalysis/${PanelApp_Date}" \
  --access-level test \
  reanalysis/reanalysis_wrapper.py \
    --conf gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json \
    --matrix gs://cpg-acute-care-main/mt/acute-care.mt \
    --pap_date "${PanelApp_Date}" \
    --ped gs://cpg-acute-care-test/reanalysis/acute_care_singleton_pedigree.json
```

---

## Structure

This program is structured using one wrapper script (`reanalysis_wrapper.py`) which defines and submits the 
worker sub-batches.

![uml](packages.png)

Optional (outputs of these scripts can be provided as files):
1. `generate_pedigree.py` - queries the Sample-Metadata API to obtain a sample pedigree for the analysis cohort
2. `process_seqr_metadata.py` - takes cohort metadata from Seqr, outputs a CPG -> Seqr ID lookup

Non-Optional
1. `query_panelapp.py` - Uses requests to obtain panel data from AU PanelApp instance
2. `hail_filter_and_classify.py` - loads variant data in Hail Query, Annotates with VEP and custom categories
3. `SLIVAR` - executed directly as commands within a Slivar container, identifies compound-het variant pairs
4. `validate_classifications.py` - takes all previous results, and runs MOI compatibility tests on all samples & variants in the cohort


---

# Module Descriptions

## Generate_pedigree.py

Utilises the Sample-Metadata API client library, which requires being executed inside a gcloud authenticated container 

Queries API for all samples in the cohort, as well as the External ID -> CPG ID mapping

Overwrites all member IDs with the CPG internal ID, and exports in PED format

**Optional argument will strip out all family structure and renders as singletons. This is to prevent the Slivar 
stage from inferring family structure whilst flagging compound hets whilst we are focused on a Singleton-only 
analysis MVP

## Process_seqr_metadata.py

A temporary measure - Seqr assigns novel IDs to families and individuals when they are inserted into the application. 
It doesn't currently expose an endpoint which can be used to map sample IDs to Seqr IDs. Once there is a viable 
endpoint to obtain this data, this mapping will not be required

We take a JSON dump of the Seqr internal sample mapping, and digest it into a simple CPG sample ID to Seqr Family 
ID map.

This remains optional throughout the analysis, but if present will allow hyperlinks from individual categorised 
variants through to the variant page for the relevant family in Seqr

## Query_panelapp.py

PanelApp is a crowd-sourced Gene-Disease knowledge-base. We are basing this analysis specifically on the 
[Mendeliome Panel](https://panelapp.agha.umccr.org/panels/137/)

This script queries for the latest version of the panel, and obtains all Green (High Confidence) genes, with 
their corresponding Mode Of Inheritance, Gene Symbol, etc. 

This data is written into the batch for downstream use as a JSON file

**Optionally takes a YYYY-MM-DD date argument. The active panel version at that date will be obtained, and the 
latest data will be augmented with differences between the latest version and this point in time (MOI change, 
or if the gene is newly Green since that date)

## Hail_filter_and_classify.py

A multi-step filtration and annotation process, using Hail Query in a PySpark cluster

Probably deserves a separate readme

## Slivar

Uses [Slivar](https://github.com/brentp/slivar) to parse the categorised VCF and search for compound-heterozygous hits

Slivar was chosen as a drop-in solution to run compound-heterozygous variant checks, respecting arbitrary family 
structures in a provided PED file. When we upgrade from Singletons to Families this can find phased comp-hets with 
no code changes required.

Implementing Comp-Hets in Hail is currently beyond me, and custom code would require corresponding development 
and testing

## Validate_classifications.py

Multi-step process for pulling in all data, then iterating through to identify variants of interest

1. Digest the Comp-Het VCF, generating a map for every sample of {Var_1: Var_2} for all variant pairs. See unit tests for output example
2. Parse the PanelApp data, and for each unique Mode of Inheritance create a filter instance. Separate Doc
3. Review each variant in turn, applying the Mode(s) of Inheritance relative to the gene the variant sits in, and the sample(s) with variant calls
4. For each sample in the dataset, record a list of all variants which passed the MOI filters
5. Pass the results to a presentation module


## MOI Tests

Flexible implementation for calculating MOI.

Includes one controlling class `MoiRunner`, and one functional class `BaseMoi` and its children. The `BaseMoi` 
derivatives each define a Mode of Inheritance.

e.g.
- DominantAutosomal - requires a variant to be exceptionally rare, and homozygotes to be absent from population 
databases. All samples with a heterozygous variant call are passed as fitting with the MOI
- XRecessive - Male samples are passed so long as the AF is below a stringent threshold, Female samples must be 
Homozygous or supported in a compound-het

The usage paradigm is:

1. Create an instance of the `MoiRunner`, passing a target MOI string and some configuration parameters
2. During the setup, the MOI string is used to determine which filters are applied to a filter list
3. The `MoiRunner` implements a `.run()` method, which takes a single variant and passes through each variant in turn
4. The result of `.run()` is a list of valid modes of inheritance, each including details of the variant(s), and samples