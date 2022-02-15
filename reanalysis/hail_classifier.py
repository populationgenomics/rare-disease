"""
- read VCF into MT
- hard-filter
- annotate
- consequence filter
- extract fields
- annotate with classes
- write as VCF
"""

from typing import Any, Dict
import json
from google.cloud import storage
import hail as hl
import click


# set some Hail constants
MISSING_STRING = hl.str('missing')
MISSING_INT = hl.int32(0)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)

# consequences we're happy to dismiss
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


def annotate_class_1(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the Class1 flag where appropriate
    :param matrix:
    :return:
    """

    return matrix


def read_json_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    this loop can be used to read config files, or data
    :param bucket_path:
    :return:
    """

    # split the full path to get the bucket and file path
    bucket = bucket_path.replace('gs://', '').split('/')[0]
    path = bucket_path.replace('gs://', '').split('/', maxsplit=1)[1]

    # create a client
    g_client = storage.Client()

    # obtain the blob of the data
    json_blob = g_client.get_bucket(bucket).get_blob(path)

    # download the blob as bytes
    return json.loads(json_blob.download_as_bytes())


def annotate_using_vep(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    hard filter some variants, and return the result
    hard filter prior to annotation to reduce runtime
    :param matrix_data:
    :return:
    """
    # hard filter for quality and abundance in the joint call
    matrix_data = matrix_data.filter_rows(matrix_data.info.AC <= 20)
    matrix_data = matrix_data.filter_rows(matrix_data.filters.length() == 0)

    # filter to biallelic loci only, and no missing variants
    matrix_data = matrix_data.filter_rows(hl.len(matrix_data.alleles) == 2)
    matrix_data = matrix_data.filter_rows(matrix_data.alleles[1] != '*')

    # throw in a repartition here (annotate even chunks in parallel)
    matrix_data = matrix_data.repartition(150, shuffle=False)

    # now run VEP 105 annotation
    vep = hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')
    return vep


def apply_row_filters(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    hard filters applied prior to the separation per-consequence
    :param matrix:
    :return:
    """
    # exac and gnomad must be below 1% or missing
    matrix = matrix.filter_rows(
        ((matrix.exac.AF < 0.01) | (hl.is_missing(matrix.exac.AF)))
        & (
            (matrix.gnomad_genomes.AF < 0.01)
            | (hl.is_missing(matrix.gnomad_genomes.AF))
        )
    )
    return matrix


def apply_consequence_filters(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    filter on per-consequence annotations
        - must be (green) genic
        - must be consequential
        - MANE transcript consequence
    :param matrix:
    :param green_genes:
    :return:
    """

    # must have a Gene ID in the green gene list
    matrix = matrix.filter_rows(
        green_genes.contains(matrix.vep.transcript_consequences.gene_id)
    )

    # require a MANE annotation
    matrix = matrix.filter_rows(
        hl.is_missing(matrix.vep.transcript_consequences.mane_select), keep=False
    )

    # discard if there are only 'useless' consequences
    matrix = matrix.filter_rows(
        hl.len(
            hl.set(matrix.vep.transcript_consequences.consequence_terms).difference(
                USELESS_CONSEQUENCES
            )
        )
        == 0,
        keep=False,
    )
    return matrix


def extract_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields from VEP and other locations
    store these in INFO
    required to be included in VCF export
    :param matrix:
    :return:
    """
    # filter the matrix table on per-consequence basis
    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            COMPOUND_CSQ=matrix.vep.transcript_consequences.gene_id
            + '|'
            + matrix.vep.transcript_consequences.transcript_id
            + '|'
            + hl.delimit(
                matrix.vep.transcript_consequences.consequence_terms, delimiter='&'
            ),
            csq=hl.delimit(
                matrix.vep.transcript_consequences.consequence_terms,
                delimiter='&',
            ),
            lof=hl.or_else(matrix.vep.transcript_consequences.lof, MISSING_STRING),
            lof_info=hl.or_else(
                matrix.vep.transcript_consequences.lof_info, MISSING_STRING
            ),
            gene_id=matrix.vep.transcript_consequences.gene_id,
            mane_transcript=matrix.vep.transcript_consequences.mane_select,
            transcript_id=matrix.vep.transcript_consequences.transcript_id,
            exon=matrix.vep.transcript_consequences.exon,
            biotype=matrix.vep.transcript_consequences.biotype,
            hgvsc=matrix.vep.transcript_consequences.hgvsc,
            hgvsp=matrix.vep.transcript_consequences.hgvsp,
            polyphen_pred=hl.or_else(
                matrix.vep.transcript_consequences.polyphen_prediction, MISSING_STRING
            ),
            polyphen_score=hl.or_else(
                matrix.vep.transcript_consequences.polyphen_score,
                MISSING_FLOAT_LO,
            ),
            sift_prediction=hl.or_else(
                matrix.vep.transcript_consequences.sift_prediction, MISSING_STRING
            ),
            sift_score=hl.or_else(
                matrix.vep.transcript_consequences.sift_score, MISSING_FLOAT_HI
            ),
            impact=matrix.vep.transcript_consequences.impact,
            exac_af=hl.or_else(matrix.exac.AF, MISSING_FLOAT_LO),
            gnomad_ex_cov=hl.or_else(matrix.gnomad_exome_coverage, MISSING_FLOAT_LO),
            gnomad_ex_af=hl.or_else(matrix.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_cov=hl.or_else(matrix.gnomad_genome_coverage, MISSING_FLOAT_LO),
            gnomad_af=hl.or_else(matrix.gnomad_genomes.AF, MISSING_FLOAT_LO),
            splice_ai_delta=hl.or_else(matrix.splice_ai.delta_score, MISSING_FLOAT_LO),
            splice_ai_csq=hl.or_else(
                matrix.splice_ai.splice_consequence, MISSING_STRING
            ).replace(' ', '_'),
            revel=hl.float64(hl.or_else(matrix.dbnsfp.REVEL_score, '0.0')),
            cadd=hl.or_else(matrix.cadd.PHRED, MISSING_FLOAT_LO),
            clinvar_sig=hl.or_else(
                matrix.clinvar.clinical_significance, MISSING_STRING
            ),
            clinvar_stars=hl.or_else(matrix.clinvar.gold_stars, MISSING_INT),
            # these next 3 are per-transcript, with ";" to delimit
            # pulling these annotations into INFO with ";" to separate
            # will break INFO parsing for most tools
            mutationtaster=hl.or_else(
                matrix.dbnsfp.MutationTaster_pred, MISSING_STRING
            ).replace(';', ','),
            fathmm=hl.or_else(matrix.dbnsfp.FATHMM_pred, MISSING_STRING).replace(
                ';', ','
            ),
            metasvm=hl.or_else(matrix.dbnsfp.MetaSVM_pred, MISSING_STRING).replace(
                ';', ','
            ),
            phast_cons=hl.float64(
                hl.or_else(matrix.dbnsfp.phastCons100way_vertebrate, '0.0')
            ),
            gerp_rs=hl.float64(hl.or_else(matrix.dbnsfp.GERP_RS, '0.0')),
            eigen_phred=hl.or_else(matrix.eigen.Eigen_phred, MISSING_FLOAT_LO),
        )
    )
    return matrix


@click.command()
@click.option('--mt', 'mt_path', help='path to the matrix table to ingest')
@click.option(
    '--ref',
    'reference',
    help='Genomic Reference to initiate Hail with',
    default='GRCh38',
)
@click.option('--pap', 'panelapp_json', help='bucket path containing panelapp JSON')
def main(mt_path: str, reference: str, panelapp_json: str):
    """

    :param mt_path:
    :param reference:
    :param panelapp_json:
    :return:
    """

    # point at the acute-care all-samples MT
    # mt_path = 'gs://cpg-acute-care-test/vep/acute-care_full_vep_105.mt'
    hl.init(default_reference=reference)

    # read MT (with annotations)
    matrix = hl.read_matrix_table(mt_path)

    matrix = annotate_using_vep(matrix_data=matrix)

    # read the parsed panelapp data from a bucket path
    panelapp = read_json_from_path(panelapp_json)

    # cast keys panel data keys (green genes) as a set of strings
    green_genes = hl.literal(set(panelapp['panel_data'].keys()))

    # explode across consequences (new row per consequence)
    matrix = matrix.explode_rows(matrix.vep.transcript_consequences)

    matrix = apply_consequence_filters(matrix, green_genes)

    # pull the annotations from 'vep' into 'INFO'
    matrix = extract_annotations(matrix)

    _matrix = annotate_class_1(matrix)

    # extract fields & and replace with placeholder if empty
    # placeholder values should be lowest possible consequence
    # e.g. most tools score 0, but for Sift 1 is least important

    # lof and lof_info are removed - always empty in current schema


# hell... maybe add some classifications while we're here...
# once the thresholds are finalised

# # export VCF, inc. Tabix index
# hl.export_vcf(
#     matrix,
#     'gs://cpg-acute-care-test/acute-care-all-vep-hail-annotated-for-slivar.vcf.bgz',
#     tabix=True,
# )
