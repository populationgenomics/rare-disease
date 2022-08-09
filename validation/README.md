# Validation

This document is notes on what has been done so far, and is not intended for `main` in its current form.

## Input

This is designed to take variant input in the form of a MatrixTable, as it is coupled with the core CPG pipeline which has MT as its main output data.

In addition to the variant data, it also sources a number of items from Metamist:

- Path to a 'Truth' VCF
- Path to a corresponding BED of confident regions
- Mapping of External:CPG IDs for the cohort

### Variant Data

A MT is provided as input, but a VCF is required for the tools to operate on. The first step is to pull the validation samples from Metamist, pull the samples from the MT, and for all common samples create a single-sample VCF. This is done with the following considerations (per sample):

- subset the MT columns to the single sample, then remove all sites where the remaining sample is HomRef
- write the remaining sites as a VCF.bgz to a `single_sample_vcf` directory, with the file containing the sample ID. These files will be persisted, and the locations recorded in Metamist for the analysis run
- `Undecided`: should all sites not indicated as PASS be removed? These are ignored by the validation software by default, so currently are being kept for debugging results.

### Metamist

Metamist is the database for the validation operations, linking each sample with the relevant data.

1. Sample IDs
    - Within metamist, samples are identified by both the external and internal IDs. Within the MT dump, the samples are referred to by their External IDs, whilst query operations on metamist objects are based on the CPG IDs. The first Metamist query is to pull a mapping of the Internal:External IDs through the `get_external_participant_id_to_internal_sample_id` method. Using this we can swap between IDs interchangeably.
2. Analysis truth
    - When running the variant data validation we need to locate the reference truth data; both the expected variant calls, and region file representing confidently called regions. In this context, the 'truth' data is treated as an Analysis result (it was not analysed by the CPG, but is an analysis product). Analysis entries have been created to store this information - see below.

#### Analysis objects

The Analysis entries in Metamist are relevant at 2 points:

1. Sourcing the data to use as a base truth
2. Recording the results of the analysis

Truth data is in the form of an expected set of variant calls, and a BED file. As well as being able to pull this data, we also need the ability to version the truth if the call-set develops over time. Metamist provides for both these requirements:

1. Analysis objects are atomic, and adding a new analysis object does not replace a previous one (though this is possible through patching); this means we can have multiple 'truth' analyses per sample at any one time
2. Analyses have a Type from a strict enumeration, so we can store these as `custom` analyses, and they will not clash with main-pipeline records
    -we could extend this enumeration to include a specific `validation` type if required
3. Analyses have an attached `Active` Boolean flag - if a truth set is replaced it would be set to inactive rather than deleted

The `meta` attribute in the Analysis objects permits storage of arbitrary data, so we can store any objects relevant to the anaylsis without creating more granular objects.

An Analysis entry has been created for each sample, with `output` indicating the analysis product (the truth VCF) and `analysis.meta` containing a pointer to the confident region BED file to use:

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
        type=AnalysisType('custom'),
        status=AnalysisStatus('completed'),
        output='gs://cpg-validation-test/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz',
        meta={
            'confident_region': 'gs://cpg-validation-test/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.bed',
        },
    ),
)
```

Note: When the same reference data is contained by multiple samples, the same analysis object can apply to all. The example here is for HG001, where we have NIST-provided and KCCG-provided sequence data, both linked to the same truth.

When processing an individual sample, we query Metamist to find the relevant truth Analysis object, then search for the VCF and BED within the Meta dictionary. Queries for these analysis objects are framed as `AnalysisQueryModel` objects, restricting to the active object(s):

```python
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_type import AnalysisType

a_query_model = AnalysisQueryModel(
    projects=['validation'],
    sample_ids=['CPG_ID'],
    type=AnalysisType('custom'),
    active=True,
)
```

Validation in the script confirms that we have exactly one `active` analysis for each sample.

Recording the results of the analysis would be done in a similar way - a new object would be created per sample, containing the relevant information:

- source MT (and/or details of the pipeline run which generated it - date/versions/other samples in callset)
- single sample VCF generated
- truth file used
- region file used
- results files generated
- top-line SENS/SPEC stats obtained

This analysis object would be of type `QC`, to differentiate the analysis inputs from analysis outputs.

## Validation Process

The validation is being carried out using [Hap.py](https://github.com/Illumina/hap.py) and [RTG's vcfeval](https://github.com/RealTimeGenomics/rtg-tools), both installed within the CPG Docker image `happy-vcfeval`. This uses Hap.py to conduct the comparison, whilst using vcfeval as the comparison engine (making use of its superior variant normalisation, in keeping with best practice advice)



## Hiccups

Truth re-processing at runtime - `vcfeval` normalises both the truth and query VCF data at runtime, and there have been a few barfs due to positional sorting. As far as I can work out, there is an issue with massive deletions - when the REF allele is LONG (e.g. hundreds of bases long) and the ALT is a single base (so a correctly represented deletion), something weird happens and when the results are re-sorted and Hap.py/VcfEval writes variants ordered by normalised position instead of POS, so BCFtools cannot sort the file. For now I've removed two variants from the Syndip truth dataset which are causing problems, e.g. the error `unsorted positions on sequence #1: 21430841 followed by 21430785` resulting from:

```chr1  21430785  .  GTTCCCAGCAGCCGCGCAGGTTTCCCCACTGGCTGCAATGGCCCTACTAAAAGCCACGTTGCATATCCGTTGTAAGCACGTGCCCTGTGCCCTGTCCCCATTCCTTATGCCCTAGGAGGCCAAGCTGGTGTCTCTAGGAGGGCCCACACGGGCACCCTGGATCCCCCAGAGAGCAGATTGGTGTGCTCAGGCCGCAGGCTGACTCAGAGCTAGGGCAGTGGGCTCTGCAGGCCACCTGGCTGGGGTTGGTGGGGGTCCTCTCTCCTGCCCCAGCTTCCACTCAGCCACCACAGTTGCCCCACCATGGGGTGGAGACGTGGGTCACCACCGGCTTGGGAGCAAGCGCCTTCTGCAGCACAGGAAGCCGAAGCTGGGGTCAGATGAGGTTGCTACCCCTGGAGGTCTGGCATAAGGGCCCCACCCTCAGGTCTCCTACACTGGCCCCATTTTACTTTGGGGTCCAAGGACAGGATGGTCAACAGGGCAGGGTGGACAGCGTGCCAGCGCCGCGCAGGGCCACCTCCCTGGGTGGATGCATCACACTAAGGAAGCGAGTGCCAAGGGGATTTAGTGGTGTGGTTCTTTCAAAGGGAGGTCAGGGTAAATGGGAATCTGCTCGGACACTCAACATGGGGGTGGGTGCACTCCTTGGAGGAGGAGGAACATGTTCAGGGGATCGTGAGGTCTTGCACAAGCCACGTGGGGCACCTTGGCTTCCCGGCAGGAGGTGGACACCCAGCCAGAGGCCTGGCTCAAGGTGACCTCACCTTCACCATGGGCTTCCTGGGTGCGCGGGCCTGAGCGCAGGTTGTTTTGTACATATTGGAATATGTGTTAACTTATGCCCCGCATCCCAACTCACACGGAAGCACGGGTCTTGTCTCAGTCTCTTCGCTGCATTTGGAAAGCAGTCTCCTCTCGGGCCAGCGCCGGGCTGAGGTGTCCAGAGGCGGCGGCAGCTGGCAGTGCCCTCAGCCCCCAAGTGTCCAGCCTGGCACTTCCCATTCAGGCCACCTGCTTTGGGTCAACAGTTCCTTTGCCAGCAGCATCTCCTAAATTGTAAGGACTCTGTCCACCCGGGGTCCTCCCAGGGCTGTGAGGACGGAAACAGGCAGGCAGTGGAGCTAACAGCTTAGTCACCAGGACCCCCAGACCTGCAAACGTCCCCTCCTGGAAGGGGAAGCCAGGAACAGCAGAACTGCCCACAAAACAAGGCTGTGAACTTTTCGGGAACTGGAACTGTTTAACTTGAACCCAGGATTGTTTAAAGCTTTATTTATTTATAGATTCTTCTTAAAAAAAAAAGTGTTGAAGAAATTTTTTGGTTATTCATACAAAAAAATGAGTTGATGATGGAAAAGCAAGTCATAATCATCTAATTGTTTTTGTCTAGGTCAAGAATGAATGTTAGCTGATAAAATAAACCCTGACAAGGAAAAAAAAAAAAAAAAAGAAATTTTACCAGACGCATATACAAAAAAGTGTGTGTGTGTGTGTGTGTGTGATCTGCACTATTTGCTTATATAAAATAACCCTGGAAGGTCACTTAAGAAACTGATAGGGCAGGCACGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATCACCTGAGGACAGGAGTTCAAGACCACCCTGGCCACCATGGCAAAACCCCTCCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGCTCACGCCTGTAATCCCAGCTACTTGGGAGGCGGAGGTGGGAGAATCGCTTGAACCCGGGAGGCAGAGGTTGCAGTGAGCCAAGATTGTGCCACTGCACTCCAGCCTGAGAGACAGAGCGAGACTCCATTTCAAAAAAAAGAAAAAGAAAAAGAAAAAAGAAACTGATAATGCTTCTCATGGTAACTTGGTGGCCAAAAAAGAAAAAGACAAAAAAAAAAAAAGAAAAAAATTGTATACATAAAAAGAAAAAGAAGCTGATAATGCTGCTTGCCTCCAAAGAAGGAAGGGGAATATCTGGGCTCAGAGGTGAGAGGGAGCCGTTTTACTGTGTACCCTTTTGAAATGTGAAACTTGTGAACATGTGACATAATCAGAAGTTAAATATAAAAATTATAAAGGCTACTTAAAATTTTGGTTTTGTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGTGCGATCTCGGCTCACTGCAACCTCTGCCTCTCTGGTTCAAGTGATTCTCCTGCCTCAGCCTCCAGAGGAGCTGGAATTACAGACGTACATCACCACACCCAGATAATTTTTGTATTTTTAGTACAGACAGGGTTTCACCATGTTGGCCAGGTTGGTCTTAAACTCCTGACCTCAGGTGATCCTCCCGCCTCAGCCTCCCAAAGTGCTGGGATTACAGGCATGAGTCACCGTGCCTGCCCTATCAGTTTCTTAAGTGACCTTCCAGGGTTATTTTATATAAGCAAATAGTGCAGATGACACACACACACACACACACACACACACACTTTTTTGTATATGCGCCTGGTAAAATTTCATCTTCACAGTTCTTGAAGATGCTTTCTTCCCCTACTAGTTATAACTGGGACGTGGCTCTGGATCAGGCCAGGAAGAACTTCTCCACTGGTTTTTCTTGGCCACAGTTTTCCACTTAAAGGATGTACCATTGTGTACTTAGCAGCAGTCCCCTTTTGATAGACATTTAGGTTGTCACCAGGATTTCTTACATGAGCAATGCTGCAGTGAACCTCCTTGTCCACAAGCCATCTCCTCTGTGCTGGGGTATATCTGTAGGATACATTCTTAGAGGAGGAGCCGGTGGGGAAGGTCCTATGCATGCATAGTTTTCGGAGCTCTCCCTCAATGCCCTCCGTTGTGGTTGTACTCACC  G  30  PASS  .  GT:AD  0|1:1,1```

File formats: VcfEval only recognises variant data with the extension `vcf.gz`, and Hail will only write a VCF file with the extension `vcf.bgz` - this means at runtime I'm running a copy to a new extension before using the files as input. Running a check now to see if this limitation is also present when running using Hap.py.
