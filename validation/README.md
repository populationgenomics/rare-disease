# Validation

## Input

Takes variant input data in the form of a MatrixTable. In addition to the variant data, it also sources per-sample
inputs from Metamist:

- Path to a 'Truth' VCFs
- Path to a corresponding BED of confident regions

### Variant Data

A MT is provided as input, but a VCF is required for the tools to operate on. The first step is to pull the validation
samples from Metamist, pull the samples from the MT, and for all common samples create a single-sample VCF. This is done
with the following considerations (per sample):

- subset the MT columns to the query single sample, then remove all sites where the remaining sample is HomRef
- write the remaining sites as a VCF.bgz to a `single_sample_vcf` directory, with the filename containing the sample ID.
    These files will be persisted, and the locations recorded in Metamist for the analysis run

### Metamist

Metamist is the database for the validation operations, linking each sample with the relevant data.

1. Sample IDs
    - Within metamist, samples are identified by both the external and internal IDs. When presenting results it will be
      beneficial to map these CPG IDs back to their External IDs for ease of reading.
2. Analysis truth
    - When running the variant data validation we need to locate the reference truth data; both the expected variant
    calls, and region file representing confidently called regions. In this context, the 'truth' data is treated as an
    Analysis result (it was not analysed by the CPG, but is an analysis product). Analysis entries have been created to
    store this information - see below.

#### Analysis objects

The Analysis entries in Metamist are relevant at 2 points:

1. Sourcing the data to use as a base truth
2. Recording the results of the analysis

Truth data is an expected set of variant calls, and a BED file. As well as being able to pull this data, we also need
the ability to version the truth if the callset develops over time. Metamist provides for these requirements:

1. Analysis objects are atomic, and adding a new analysis object does not replace a previous one (though this is
    possible through patching); this means we can have multiple 'truth' analyses per sample at any one time
2. Analyses have an attached `Active` Boolean flag - if a truth set is replaced it would be set to inactive rather than
    deleted; this active/inactive state can be included in the metamist queries.
3. Analyses have a Type from a strict enumeration, so we store these as `qc` analyses, with a sub-type added in the
    entry `meta`. We can phrase queries on these objects using the meta values, allowing us flexibility in how we query
    for and store result sets.

An Analysis entry has been created for each sample, with `output` indicating the analysis product (the truth VCF) and
`analysis.meta` containing a pointer to the confident region BED file to use. If this is replaced/updated for any
sample, the existing analysis entry should be set as inactive, then the new one added (no deletions):

```python
from sample_metadata.apis import AnalysisApi
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_status import AnalysisStatus

anal_api = AnalysisApi()
anal_api.create_new_analysis(
    project='validation',
    analysis_model=AnalysisModel(
        sample_ids=['CPG_ID', 'CPG_ID2'],
        type=AnalysisType('qc'),
        status=AnalysisStatus('completed'),
        output='gs://cpg-validation-test/HG001/HG001_GRCh38_benchmark.vcf.gz',
        meta={
            'type': 'validation',
            'confident_region': 'gs://cpg-validation-test/HG001/HG001_GRCh38_benchmark.bed',
        },
    ),
)
```

Note: When the same reference data is contained by multiple samples, the same analysis object can apply to all. An
example is NA12878, where we have NIST-provided and KCCG-provided sequence data, both linked to the same truth.

When processing an individual sample, we query Metamist to find the relevant truth Analysis object, then search for the
VCF and BED within object data. Queries for these analyses are framed as `AnalysisQueryModel` objects:

```python
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_type import AnalysisType

a_query_model = AnalysisQueryModel(
    projects=['validation'],
    sample_ids=['CPG_ID'],
    type=AnalysisType('qc'),
    active=True,
)
```

Validation in the script confirms that we have exactly one `active` analysis for each sample - if there are multiple
active sets of truth data, the stage will fail with an informative message. Recording the results of the analysis is
done in a similar way - a new object is created per sample, containing relevant information:

- source MT (and/or details of the pipeline run which generated it - date/versions/other samples in callset)
- location of the dated comparison results folder
- single sample VCF generated
- truth VCF used
- truth BED used
- results files generated
- top-line SENS/SPEC stats obtained (parsed from the summary CSV)

## Validation Process

The validation is being carried out using [Hap.py](https://github.com/Illumina/hap.py) and
[RTG's vcfeval](https://github.com/RealTimeGenomics/rtg-tools), both installed within the CPG Docker image `hap.py`.
Hap.py is used to conduct the varaint comparison, with vcfeval as the comparison engine (making use of its superior
variant normalisation, in keeping with best practice advice)
