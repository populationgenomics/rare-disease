"""
script to be used in vertex notebook for generating
a slivar-friendly file for use in queries

Even though we are exploding the rows, which should inflate
the overall VCF size, we are retaining only a subset of transcripts

effectively, this is the same as having only annotated against MANE
transcripts in the first place, then using split-vep to bring all
those annotations into the INFO field of each row

purely in terms of file size, this is WAY smaller than retaining
all the redundant transcript consequences, and can be done prior to
writing any data out to disk

---- work out what to do with this. We can run an end-to-end
- read VCF into MT
- hard-filter
- annotate
- filter
- extract fields
- write as VCF

OR
this can be a standalone process, using an annotated MT

IMO this is probably useful with a checkpoint, e.g. read INPUT
(VCF or MT)
annotate, save in full as MT, extract & filter, save as VCF
"""


import hail as hl


# point at the acute-care all-samples MT
MT_PATH = 'gs://cpg-acute-care-test/vep/acute-care_full_vep_105.mt'
hl.init(default_reference='GRCh38')

# read MT (with annotations)
matrix = hl.read_matrix_table(MT_PATH)

# do some filtering
# either one is below the 1% threshold, or both are missing
matrix = matrix.filter_rows(
    ((matrix.exac.AF <= 0.01) | (matrix.gnomad_genomes.AF <= 0.01))
    | ((hl.is_missing(matrix.exac.AF)) & (hl.is_missing(matrix.gnomad_genomes.AF)))
)

# max 20 times called within this  169 sample joint call. Still V. high
matrix = matrix.filter_rows(matrix.info.AC <= 20)

# strict filter on FILTER
matrix = matrix.filter_rows(matrix.filters.length() == 0)

# remove some row-annotation redundancy (this could affect the object size)
matrix = matrix.drop('mainTranscript', 'sortedTranscriptConsequences')

# explode across consequences (new row per consequence)
matrix = matrix.explode_rows(matrix.vep.transcript_consequences)

USELESS_CONSEQUENCES = hl.literal(
    {
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "downstream_gene_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_exon_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "mature_miRNA_variant",
    }
)

# filter some per-cons (e.g. must be genic, must be consequential, MANE transcript)
# MANE filtering is on MANE_SELECT, MANE_PLUS_CLINICAL not yet in the VEP-hail schema
matrix = matrix.filter_rows(
    (hl.is_missing(matrix.vep.transcript_consequences.gene_id))
    | (
        hl.len(
            hl.set(matrix.vep.transcript_consequences.consequence_terms).difference(
                USELESS_CONSEQUENCES
            )
        )
        == 0
    )
    | (hl.is_missing(matrix.vep.transcript_consequences.mane_select)),
    keep=False,
)

# extract fields we're interested in, and replace crucial values with
# placeholder values if they are currently empty
# placeholder values should be lowest possible consequence
# e.g. most tools score 0, but for Sift 1 is least important
matrix = matrix.annotate_rows(
    info=matrix.info.annotate(
        vep_csq=hl.delimit(
            matrix.vep.transcript_consequences.consequence_terms,
            delimiter='&',
        ),
        vep_gene_ids=matrix.vep.transcript_consequences.gene_id,
        vep_mane_select=matrix.vep.transcript_consequences.mane_select,
        vep_transcript_id=matrix.vep.transcript_consequences.transcript_id,
        vep_canonical=matrix.vep.transcript_consequences.canonical,
        vep_exon=matrix.vep.transcript_consequences.exon,
        vep_biotype=matrix.vep.transcript_consequences.biotype,
        vep_hgvsc=matrix.vep.transcript_consequences.hgvsc,
        vep_hgvsp=matrix.vep.transcript_consequences.hgvsp,
        vep_lof=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.lof),
            hl.str('missing'),
            matrix.vep.transcript_consequences.lof,
        ),
        vep_lof_info=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.lof_info),
            hl.str('missing'),
            matrix.vep.transcript_consequences.lof_info,
        ),
        vep_polyphen_prediction=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.polyphen_prediction),
            hl.str('missing'),
            matrix.vep.transcript_consequences.polyphen_prediction,
        ),
        # for polyphen 1 is deleterious, 0 is tolerated
        vep_polyphen_score=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.polyphen_score),
            hl.float64(0.0),
            matrix.vep.transcript_consequences.polyphen_score,
        ),
        vep_sift_prediction=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.sift_prediction),
            hl.str('missing'),
            matrix.vep.transcript_consequences.sift_prediction,
        ),
        # for SIFT, high 1 is tolerated, 0 is damaging
        vep_sift_score=hl.if_else(
            hl.is_missing(matrix.vep.transcript_consequences.sift_score),
            hl.float64(1.0),
            matrix.vep.transcript_consequences.sift_score,
        ),
        vep_impact=matrix.vep.transcript_consequences.impact,
        exac_af=hl.if_else(
            hl.is_missing(matrix.exac.AF), hl.float64(0.0), matrix.exac.AF
        ),
        gnomad_ex_cov=hl.if_else(
            hl.is_missing(matrix.gnomad_exome_coverage),
            hl.float64(0.0),
            matrix.gnomad_exome_coverage,
        ),
        gnomad_ex_af=hl.if_else(
            hl.is_missing(matrix.gnomad_exomes.AF),
            hl.float64(0.0),
            matrix.gnomad_exomes.AF,
        ),
        gnomad_cov=hl.if_else(
            hl.is_missing(matrix.gnomad_genome_coverage),
            hl.float64(0.0),
            matrix.gnomad_genome_coverage,
        ),
        gnomad_af=hl.if_else(
            hl.is_missing(matrix.gnomad_genomes.AF),
            hl.float64(0.0),
            matrix.gnomad_genomes.AF,
        ),
        splice_ai_delta=hl.if_else(
            hl.is_missing(matrix.splice_ai.delta_score),
            hl.float64(0.0),
            matrix.splice_ai.delta_score,
        ),
        splice_ai_csq=hl.if_else(
            hl.is_missing(matrix.splice_ai.splice_consequence),
            hl.str('missing'),
            matrix.splice_ai.splice_consequence.replace(' ', '_'),
        ),
        revel=hl.if_else(
            hl.is_missing(matrix.dbnsfp.REVEL_score),
            hl.str('missing'),
            matrix.dbnsfp.REVEL_score,
        ),
        cadd=hl.if_else(
            hl.is_missing(matrix.cadd.PHRED), hl.float64(0.0), matrix.cadd.PHRED
        ),
        clinvar_sig=hl.if_else(
            hl.is_missing(matrix.clinvar.clinical_significance),
            hl.str('missing'),
            matrix.clinvar.clinical_significance,
        ),
        clinvar_stars=hl.if_else(
            hl.is_missing(matrix.clinvar.gold_stars),
            hl.int32(0),
            matrix.clinvar.gold_stars,
        ),
    )
)

# hell... maybe add some classifications while we're here...

# # export VCF, inc. Tabix index
# hl.export_vcf(
#     matrix,
#     'gs://cpg-acute-care-test/acute-care-all-vep-hail-annotated-for-slivar.vcf.bgz',
#     tabix=True,
# )
