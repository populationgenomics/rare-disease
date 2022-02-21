# Reanalysis MVP

---

Multi-step Rare Disease variant prioritisation Framework

1. Queries PanelApp for the latest Mendeliome data
2. Loads variant data into a Hail MatrixTable
    * Hard-Filters variants on QC attributes
    * Annotates with VEP 105
    * Filters variants out using consequnce annotations
    * Applies classifications approximating the Miro specification
    * Filters to classified variants only, and exports to VCF
3. Runs Slivar to obtain all compound-het variants
4. Runs analysis module to create final classifications
    * Import Compound-Het VCF, and record variant pairs
    * Read classified VCF, and validate classifications (see below)
5. Report the classified variants in an HTML format

---

### Hail Stage

Hard filters for variant removal:

* Any calls failing a VCF `FILTER`
* Any non-normalised variants (> 2 alleles)
* Any missing calls (ALT = '*')
* In a large enough joint call, remove variants present in over 10% of samples as likely noise

Broad consequence-based filtering:

* Variants not in a PanelApp 'Green' gene
* Variants are Rare in Gnomad and Exac (default threshold is 1%)
* Consequences rated as `LOW` in VEP
* Clinvar Benign with at least one `star`
* If the variant has a MANE transcript consequence, filter to only retain that csq

Classification steps:

```
Note: Classes are currently non-exclusive. 
i.e. each Variant can be scored against multiple categories
An analysis step could involve showing each variant against 
only its most severe class; 1 > 2 > 3 > 4
```

Class 1: 
* Rare
* Clinvar Pathogenic/Likely Pathogenic
* At least 1 Clinvar star 

Class 2: 
* Clinvar Pathogenic/Likely Pathogenic (any stars), or 
* VEP High consequence variant annotation, or 
* Missense with pathogenic consensus in CADD & Revel
* Requires validation that MOI or gene Green status is new in PanelApp (relative to a given data, given as a parameter)

Class 3: 
* VEP High Impact consequence
* Supported by Loftee, or 
* Clinvar annotation Pathogenic/Likely Pathogenic (any stars) 
* *pending refinement, Loftee data currently not populated in the Hail data model

Class 4: 
* Rare, and either 
* CADD and Revel consensus damaging prediction, or
* Predicted damaging consensus across [MutationTaster, Sift, Polyphen, Gerp, Eigen]

NB. class assignments are done independently of MOI and genotype, as those filters are applied outside of Hail

---

### Compound Het checks

To assist the downstream MOI check, the VCF containing all variants of interest is analysed using Slivar, identifying all variant pairs forming an intra-genic compound-het. 
This 3rd party tool is preferred to creating a custom solution; for singleton analysis the check is trivial, but Slivar has already been built to recognise phased inheritance in arbitrary family structures - this is a plug-in we can easily use from the start, and transition from singleton to familial analysis without skipping a beat

---

### Class & MOI analysis

This component brings in the PanelApp, Classification VCF, and Compound-Het VCF

1. Digest the compound-het VCF, and store as a lookup of `Sample -> Variant -> Pair`
2. Digest the PanelApp data, storing a representation of the MOI to apply to each gene