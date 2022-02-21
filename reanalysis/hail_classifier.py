"""
- read VCF into MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC,
- annotate
- consequence filter
- extract fields
- annotate with classes
- write as VCF
"""

from typing import Any, Dict
import json
import logging
import sys

import click
from google.cloud import storage
import hail as hl


# set some Hail constants
MISSING_STRING = hl.str('missing')
MISSING_INT = hl.int32(0)
ONE_INT = hl.int32(1)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)


def annotate_mane(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    annotates into each row whether there was a MANE transcript present
    (only one of the transcripts has to have a MANE equivalent)
    :param matrix:
    :return:
    """

    return matrix.annotate_rows(
        mane_tx_present=hl.if_else(
            hl.delimit(
                hl.flatten(matrix.vep.transcript_consequences.consequence_terms),
                delimiter='|',
            ).contains('NM'),
            ONE_INT,
            MISSING_INT,
        )
    )


def annotate_class_1(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    applies the Class1 flag where appropriate
    semi-rare in Gnomad
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
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
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
    rare in Gnomad
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
                    (
                        hl.len(
                            hl.set(matrix.info.csq).intersection(critical_consequences)
                        )
                        > 0
                    )
                    | (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                )
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def annotate_class_2(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Provisional class, requires a panelapp check (new)
    rare in Gnomad, and
    - Clinvar, or
    - Critical protein consequence, or
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
                (
                    (
                        hl.len(
                            hl.set(matrix.info.csq).intersection(critical_consequences)
                        )
                        > 0
                    )
                    | (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                    | (
                        (matrix.info.cadd > config.get('cadd'))
                        & (matrix.info.revel > config.get('revel'))
                    )
                )
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def annotate_class_4(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Filler Class, based on in silico annotations
    CADD & REVEL above threshold (switched to consensus), or
    Massive cross-tool consensus
    rare in Gnomad
    :param matrix:
    :param config:
    :return:
    """

    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            Class4=hl.if_else(
                (
                    (
                        (matrix.info.cadd > config.get('cadd'))
                        & (matrix.info.revel > config.get('revel'))
                    )
                    | (
                        (matrix.info.sift_score < config.get('sift'))
                        & (matrix.info.polyphen_score >= config.get('polyphen'))
                        & (
                            (matrix.info.mutationtaster.contains("D"))
                            | (matrix.info.mutationtaster == "missing")
                        )
                        & (matrix.info.gerp_rs >= config.get('gerp'))
                        & (matrix.info.eigen_phred > config.get('eigen'))
                    )
                )
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    return matrix


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a GCP bucket path to a JSON file, read into an object
    this loop can read config files, or data
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

    # the download_as_bytes method isn't available; but this returns bytes?
    json_content = json_blob.download_as_string()

    return json.loads(json_content)


def hard_filter_before_annotation(
    matrix_data: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """

    :param matrix_data:
    :param config:
    :return:
    """

    # count the samples in the VCF
    # if we reach the sample threshold, filter on AC
    if matrix_data.count_cols() >= config.get('min_samples_to_ac_filter'):
        matrix_data = matrix_data.filter_rows(
            matrix_data.info.AC
            <= matrix_data.info.AN // config.get('ac_filter_percentage')
        )

    # filter to cohort samples
    # filter to variants in THESE samples

    # hard filter for quality
    # assumption here that data is well normalised in pipeline
    matrix_data = matrix_data.filter_rows(
        (
            matrix_data.filters.length() == 0
        )  # clarify with FILTER='PASS' & GATK tranches
        & (hl.len(matrix_data.alleles) == 2)
        & (matrix_data.alleles[1] != '*')
    )

    # throw in a repartition here (annotate even chunks in parallel)
    matrix_data = matrix_data.repartition(150, shuffle=False)

    return matrix_data


def annotate_using_vep(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    runs VEP annotation on the MT
    :param matrix_data:
    :return:
    """

    # now run VEP 105 annotation
    vep = hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')
    return vep


def apply_row_filters(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    variant filters applied prior to the per-consequence split
    where the population frequencies at?
    :param matrix:
    :param config:
    :return:
    """
    # exac and gnomad must be below threshold or missing
    matrix = matrix.filter_rows(
        (
            (matrix.exac.AF < config.get('exac_semi_rare'))
            | (hl.is_missing(matrix.exac.AF))
        )
        & (
            (matrix.gnomad_genomes.AF < config.get('gnomad_semi_rare'))
            | (hl.is_missing(matrix.gnomad_genomes.AF))
        )
    )

    # remove all clinvar benign, any level of support
    # this isn't working at the moment...
    matrix = matrix.filter_rows(
        (matrix.clinvar.clinical_significance.upper().contains('BENIGN'))
        & (matrix.clinvar.gold_stars > 0),
        keep=False,
    )
    return matrix


def apply_consequence_filters(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    filter on per-consequence annotations
        - must be (green) genic
        - must be consequential
        - MANE transcript consequence*

    Removing all consequence annotations without Mane might be excessive here...
    Allow for retention where the row has other attributes
    :param matrix:
    :param green_genes:
    :param config:
    :return: input matrix filtered by consequences
    """

    # identify consequences to discard from the config
    useless_csq = hl.set(config.get('useless_csq'))

    # must have a Gene ID in the green gene list
    matrix = matrix.filter_rows(
        green_genes.contains(matrix.vep.transcript_consequences.gene_id)
    )

    # filter out non-Mane transcript if the variant has at least one MANE annotation
    # require high impact consequences
    # enough csq impact is 'at least 1 csq outside useless set'
    matrix = matrix.filter_rows(
        (
            (hl.is_missing(matrix.vep.transcript_consequences.mane_select))
            & (matrix.mane_tx_present == 1)
        )
        | (
            hl.len(
                hl.set(matrix.vep.transcript_consequences.consequence_terms).difference(
                    useless_csq
                )
            )
            == 0
        ),
        keep=False,
    )

    # expecting a minority remaining - repartition
    logging.info('Repartition to 40 fragments following CSQ filter')
    matrix = matrix.repartition(40, shuffle=False)

    return matrix


def extract_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields from VEP and other locations
    store these in INFO (required to be included in VCF export)

    replace with placeholder if empty
    placeholder values should be least consequential
    e.g. most tools score 0, but for Sift 1 is least important

    Adds in a compound CSQ field, useful with Slivar
    :param matrix:
    :return: input matrix with annotations pulled into INFO
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
            polyphen_prediction=hl.or_else(
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
    :return: input matrix, minus rows without Classes applied
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
    """
    hl.export_vcf(
        matrix,
        output_path,
        tabix=True,
    )


@click.command()
@click.option('--mt', 'mt_path', help='path to the matrix table to ingest')
@click.option('--pap', 'panelapp_path', help='bucket path containing panelapp JSON')
@click.option('--config', 'config_path', help='path to a config dict')
@click.option('--output', 'out_vcf', help='VCF path to export results')
def main(mt_path: str, panelapp_path: str, config_path: str, out_vcf: str):
    """
    Read the MT from disk, do filtering and class annotation
    Export as a VCF

    :param mt_path: path to the MT dump
    :param panelapp_path: path to the panelapp data dump
    :param config_path: path to the config json
    :param out_vcf: path to write the VCF out to
    """

    logging.info('Reading config dict from "%s"', config_path)
    # get the run configuration JSON
    config_dict = read_json_dict_from_path(config_path)

    logging.info('Reading PanelApp data from "%s"', panelapp_path)
    # read the parsed panelapp data from a bucket path
    panelapp = read_json_dict_from_path(panelapp_path)

    # cast panel data keys (green genes) as a set(str)
    green_genes = set(panelapp.keys())
    green_gene_set_expression = hl.literal(green_genes)
    logging.info('Extracted %d green genes', len(green_genes))

    logging.info(
        'Starting Hail with reference genome "%s"', config_dict.get('ref_genome')
    )
    # initiate Hail with the specified reference
    hl.init(default_reference=config_dict.get('ref_genome'), quiet=True)

    logging.info('Loading MT from "%s"', mt_path)
    # load MT in
    matrix = hl.read_matrix_table(mt_path)

    # hard filter entries in the MT prior to annotation
    logging.info('Hard filtering variants')
    matrix = hard_filter_before_annotation(matrix_data=matrix, config=config_dict)

    # re-annotate using VEP
    logging.info('Annotating variants')
    matrix = annotate_using_vep(matrix_data=matrix)

    # filter on consequence-independent row annotations
    logging.info('Filtering Variant rows')
    matrix = apply_row_filters(matrix=matrix, config=config_dict)

    # annotate the MANE transcript property on each variant
    # done prior to exploding on consequence
    matrix = annotate_mane(matrix)

    # explode consequences (new row per csq)
    logging.info('Splitting variant rows by consequence')
    matrix = matrix.explode_rows(matrix.vep.transcript_consequences)

    # pull the annotations from 'vep' into 'info'
    logging.info('Pulling VEP annotations into INFO field')
    matrix = extract_annotations(matrix)

    # hard filter for relevant consequences
    logging.info('Hard filter rows on consequence')
    matrix = apply_consequence_filters(matrix, green_gene_set_expression, config_dict)

    # add Classes to the MT
    logging.info('Applying classes to variant consequences')
    matrix = annotate_class_1(matrix, config_dict)
    matrix = annotate_class_3(matrix, config_dict)
    matrix = annotate_class_2(matrix, config_dict)
    matrix = annotate_class_4(matrix, config_dict)

    # possibly add a background class here for interesting, but only
    # good enough to be a second hit. C4 is this for now

    # Need a new section here - if we split on consequence
    # to expose more rows, we need a way to reduce that number
    # of rows back down. Hard filtering on MANE is a bit...
    # and NOT doing that means that the less stringent categories
    # can have 10 repetitions of the same variant

    # maybe if a gene has a locus has a mane transcript, crush it
    # if it doesn't, then don't

    # filter to class-annotated only prior to export
    logging.info('Filter variants to leave only classified')
    matrix = filter_to_classified(matrix)

    # write the results to a VCF path
    logging.info('Write variants out to "%s"', out_vcf)
    write_matrix_to_vcf(matrix=matrix, output_path=out_vcf)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
