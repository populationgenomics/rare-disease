# Genomic Region Stratification

[This folder](gs://cpg-validation-test/GRCh38_regions) contains the BED files used to further stratify the variant calling validation process.

## Data Source

These BED files are derived from the GIAB regions, which are defined/created in
[this repository](https://github.com/genome-in-a-bottle/genome-stratifications/), and downloadable
in full from [the corresponding FTP site](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/).
The [Hap.py repository](https://github.com/Illumina/hap.py) contains a description of how these
extended regions can be integrated into a validation run
[here](https://github.com/Illumina/hap.py/blob/master/doc/happy.md#stratification-via-bed-regions).

The genome-stratifications repository provides BED regions divided by genome build, then by category.
Within these categories there is a lot of granularity. Here we are only using a subset of these in our validations.

## Files

1. `definition.tsv`, containing a short description and relative path to a BED file
2. `*.bed.gz`, compressed files containing the region definitions for each category

The `TSV` file is used by `hap.py` to locate region files, and annotate them in the output file with the short description. The files currently in use are:

- `refseq_cds`: RefSeq Coding regions only [see here](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/FunctionalRegions/GRCh38-FunctionalRegions-README.md#data--file-overview)
- `nonunique`: Mappability, regions of high homology over 250bp [see here](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/mappability/GRCh38-mappability-README.md#data--file-overview)
- `gnomAD inbreeding`: [gnomAD inbreedingcoeff variants](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/OtherDifficult/GRCh38-OtherDifficult-README.md#data--file-overview)
- `MHC`: [chr6 MHC complex regions](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/OtherDifficult/GRCh38-OtherDifficult-README.md#data--file-overview)
- `all_difficult`: Union of all complex regions (e.g. tandem repeats, all homopolymers >6bp, all imperfect homopolymers >10bp, all difficult to map regions, all segmental duplications, GC <25% or >65%, "Bad Promoters", chrX/Y XTR and ampliconic, satellites) [see](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/union/GRCh38-union-README.md#file-descriptions)
- `not_in_all_difficult`: Complement, whole genome other than `all_difficult` regions

## Usage

In practice, these files are used in the following way:

- when running the validation workflow provide the path to the stratification directory of choice (e.g. `gs://cpg-validation-test/GRCh38_regions`)
- when validation stages are being created, the script checks for a `definition.tsv` file, and copies all files into the execution container
  - N.B. it is important that all the files are copied simultaneously, otherwise Hail will copy all files into different temp folders
- within the container, the hap.py script has the argument `--stratification /location/of/definition.tsv` added

## Results

In the Hap.py output files there will be one ending with `extended.csv`, with entries for the whole (confident) ROI,
then each separate stratified location. These results will be logged into metamist for the sample(s).

I can't see where this is explicitly stated in the Hap.py documentation, but testing shows that the stratification
regions are subsets within the confident regions. i.e. we don't require separate stratification folders for each capture,
it is sufficient to update the confident regions BED, and all stratifications will be relative to that. In the hap.py
codebase there are multiple layers in multiple languages, and as far as I can see the specific documentation that all
stratified regions are subsets within the confident region is
[this](https://github.com/Illumina/hap.py/blob/master/src/c%2B%2B/include/QuantifyRegions.hh#L98)
