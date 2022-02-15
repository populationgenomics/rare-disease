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
ONE_INT = hl.int32(1)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)


def annotate_class_1(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    applies the Class1 flag where appropriate
    rare (< 0.005 in Gnomad)
    at least one Clinvar star
    either Pathogenic or Likely_pathogenic in Clinvar
    Assign 1 or 0, depending on presence
    :param matrix:
    :param config:
    :return:
    """

    clinvar_pathogenic_terms = hl.set(config.get('clinvar_path'))

    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            Class1=hl.if_else(
                (matrix.info.clinvar_stars > 0)
                & (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                & (matrix.info.gnomad_af < config.get('gnomad')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def annotate_class_3(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    applies the Class3 flag where appropriate
    at least one important consequence
    rare (< 0.005 in Gnomad)
    either predicted NMD (Loftee not in the data yet) or
    any star Pathogenic or Likely_pathogenic in Clinvar
    :param matrix:
    :param config:
    :return:
    """
    critical_consequences = hl.set(config.get('critical_csq'))
    clinvar_pathogenic_terms = hl.set(config.get('clinvar_path'))
    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            Class3=hl.if_else(
                (
                    hl.len(hl.set(matrix.info.csq).intersection(critical_consequences))
                    > 0
                )
                & (
                    (matrix.info.clinvar_stars > 0)
                    & (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                )
                & (matrix.info.gnomad_af < config.get('gnomad')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def annotate_class_2(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Provisional class, requires a panelapp check (new)
    Basically for now this is
    rare (< 0.005 in Gnomad), and
    - Class 3, or
    - Clinvar, or
    - High in silico consequence
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))
    clinvar_pathogenic_terms = hl.set(config.get('clinvar_path'))

    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            Class2=hl.if_else(
                (matrix.info.Class3 == 1)
                | (
                    hl.len(hl.set(matrix.info.csq).intersection(critical_consequences))
                    > 0
                )
                & (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                & (matrix.info.gnomad_af < config.get('gnomad')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def annotate_class_4(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Filler Class, based on in silico annotations
    CADD/REVEL above threshold, or
    Massive cross-tool consensus
    rare (< 0.005 in Gnomad)
    :param matrix:
    :param config:
    :return:
    """

    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            Class4=hl.if_else(
                (
                    (matrix.info.cadd > config.get('cadd'))
                    | (matrix.info.revel > config.get('revel'))
                )
                | (
                    (matrix.info.vep_sift_prediction == 'deleterious')
                    & (matrix.info.vep_polyphen_prediction == 'probably_damaging')
                    & (matrix.info.vep_sift_score < config.get('sift'))
                    & (matrix.info.vep_polyphen_score >= config.get('polyphen'))
                    & (
                        (matrix.info.mutationtaster.includes("D"))
                        | (matrix.info.mutationtaster == "missing")
                    )
                    & (matrix.info.gerp_rs >= config.get('gerp'))
                    & (matrix.info.eigen_phred > config.get('eigen'))
                )
                & (matrix.info.gnomad_af < config.get('gnomad')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

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
    matrix: hl.MatrixTable, green_genes: hl.SetExpression, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    filter on per-consequence annotations
        - must be (green) genic
        - must be consequential
        - MANE transcript consequence
    :param matrix:
    :param green_genes:
    :param config:
    :return:
    """
    useless_csq = hl.set(config.get('useless_csq'))

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
                useless_csq
            )
        )
        == 0,
        keep=False,
    )
    return matrix


def extract_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields from VEP and other locations
    store these in INFO (required to be included in VCF export)

    replace with placeholder if empty
    placeholder values should be least consequential
    e.g. most tools score 0, but for Sift 1 is least important
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
            csq=hl.set(matrix.vep.transcript_consequences.consequence_terms),
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


def filter_to_classified(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filters to only rows which have an associated class
    :param matrix:
    :return:
    """
    matrix = matrix.filter_rows(
        (matrix.info.Class1 == 1)
        | (matrix.info.Class2 == 1)
        | (matrix.info.Class3 == 1)
        | (matrix.info.Class4 == 1)
    )
    return matrix


def write_matrix_to_vcf(matrix: hl.MatrixTable, output_path: str):
    """

    :param matrix:
    :param output_path: where to write the VCF
    :return:
    """
    hl.export_vcf(
        matrix,
        output_path,
        tabix=True,
    )


@click.command()
@click.option('--mt', 'mt_path', help='path to the matrix table to ingest')
@click.option('--pap', 'panelapp_json', help='bucket path containing panelapp JSON')
@click.option('--config', 'config', help='path to a config dict')
@click.option('--output', 'out_vcf', help='VCF path to export results')
def main(mt_path: str, panelapp_json: str, config: str, output: str):
    """
    Read the MT from disk, do filtering and class annotation
    Export as a VCF

    :param mt_path: path to the MT dump
    :param panelapp_json: path to the panelapp data dump
    :param config: path to the config json
    :param output: path to write the VCF out to
    :return:
    """

    # get the run configuration JSON
    config_dict = read_json_from_path(config)

    # initiate Hail with the specified reference
    hl.init(default_reference=config_dict.get('ref_genome'))

    # load MT in
    matrix = hl.read_matrix_table(mt_path)

    # re-annotate using VEP
    matrix = annotate_using_vep(matrix_data=matrix)

    # read the parsed panelapp data from a bucket path
    panelapp = read_json_from_path(panelapp_json)

    # cast panel data keys (green genes) as a set(str)
    green_genes = hl.literal(set(panelapp['panel_data'].keys()))

    # explode consequences (new row per csq)
    matrix = matrix.explode_rows(matrix.vep.transcript_consequences)

    # hard filter for relevant consequences
    matrix = apply_consequence_filters(matrix, green_genes, config_dict)

    # pull the annotations from 'vep' into 'info'
    matrix = extract_annotations(matrix)

    # add Classes to the MT
    matrix = annotate_class_1(matrix, config_dict)
    matrix = annotate_class_3(matrix, config_dict)
    matrix = annotate_class_2(matrix, config_dict)
    matrix = annotate_class_4(matrix, config_dict)

    # filter to class-annotated only prior to export
    matrix = filter_to_classified(matrix)

    # write the results to a VCF path
    write_matrix_to_vcf(matrix=matrix, output_path=output)
