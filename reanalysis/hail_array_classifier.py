"""
Similar to hail_classifier.py script, but no exploding
- read VCF into MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- annotate
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- consequence filter
- remove all rows with no interesting GREEN consequences
- extract vep data into CSQ string(s)
- annotate with classes
- remove all un-classified variants
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

    return matrix.annotate_rows(
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


def annotate_class_2(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Provisional class, requires a panelapp check
    - Rare in Gnomad, and
    - Clinvar, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))
    clinvar_pathogenic_terms = hl.set(config.get('clinvar_path'))

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class2=hl.if_else(
                (
                    matrix.vep.transcript_consequences.any(
                        lambda x: hl.len(
                            critical_consequences.intersection(
                                hl.set(x.consequence_terms)
                            )
                        )
                        > 0
                    )
                )
                | (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                | (
                    (matrix.info.cadd > config.get('cadd'))
                    & (matrix.info.revel > config.get('revel'))
                )
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_class_3(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    applies the Class3 flag where appropriate
    - Critical protein consequence on at least one transcript
    - rare in Gnomad
    - either predicted NMD (Loftee not in the data yet) or
    - any star Pathogenic or Likely_pathogenic in Clinvar

    currently this class is a bit baggy, as LOF confirmation doesn't exist in Hail-VEP
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))
    clinvar_pathogenic_terms = hl.set(config.get('clinvar_path'))

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class3=hl.if_else(
                (
                    matrix.vep.transcript_consequences.any(
                        lambda x: hl.len(
                            critical_consequences.intersection(
                                hl.set(x.consequence_terms)
                            )
                        )
                        > 0
                    )
                )
                & (clinvar_pathogenic_terms.contains(matrix.info.clinvar_sig))
                & (matrix.info.gnomad_af < config.get('gnomad_semi_rare')),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_class_4(matrix: hl.MatrixTable, config: Dict[str, Any]) -> hl.MatrixTable:
    """
    Class based on in silico annotations
    - rare in Gnomad, and
    - CADD & REVEL above threshold (switched to consensus), or
    - Massive cross-tool consensus
    :param matrix:
    :param config:
    :return:
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class4=hl.if_else(
                (matrix.info.gnomad_af < config.get('gnomad_semi_rare'))
                & (
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
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


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
    return json.loads(json_blob.download_as_string())


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

    # hard filter for quality
    # assumption here that data is well normalised in pipeline
    matrix_data = matrix_data.filter_rows(
        (
            matrix_data.filters.length() == 0
        )  # clarify with FILTER='PASS' & GATK SNP tranches
        & (hl.len(matrix_data.alleles) == 2)
        & (matrix_data.alleles[1] != '*')
    )

    # throw in a repartition here (annotate even chunks in parallel)
    return matrix_data.repartition(150, shuffle=False)


def annotate_using_vep(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    runs VEP annotation on the MT
    :param matrix_data:
    :return:
    """

    # now run VEP 105 annotation
    return hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')


def filter_mt_rows(
    matrix: hl.MatrixTable, config: Dict[str, Any], green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - remove common variants
    - reduce the 'geneIds' set to contain only green genes
    - reduce the per-row transcript consequences to those specific to the geneIds
    - reduce the rows to ones where there are remaining tx consequences
    :param matrix:
    :param config:
    :param green_genes:
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

    # remove all clinvar benign, decent level of support
    benign = hl.str('benign')
    matrix = matrix.filter_rows(
        (matrix.clinvar.clinical_significance.lower().contains(benign))
        & (matrix.clinvar.gold_stars > 0),
        keep=False,
    )

    # remove any rows with no genic annotation at all
    # do this here as set intersections with Missing will fail
    matrix = matrix.filter_rows(hl.is_missing(matrix.geneIds), keep=False)

    # replace the default list of green IDs with a reduced set
    matrix = matrix.annotate_rows(geneIds=green_genes.intersection(matrix.geneIds))

    # split to form a separate row for each green gene
    # this transforms the 'geneIds' field from a set to a string
    matrix = matrix.explode_rows(matrix.geneIds)

    # identify consequences to discard from the config
    useless_csq = hl.set(config.get('useless_csq'))

    # reduce consequences to overlap with per-variant green geneIDs (pre-filtered)
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: (matrix.geneIds == x.gene_id)
                & (hl.len(hl.set(x.consequence_terms).difference(useless_csq)) != 0)
            )
        )
    )

    # filter out all rows with no remaining consequences
    matrix = matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) > 0)

    logging.info('Repartition to 50 fragments following Gene ID filter')
    matrix = matrix.repartition(50, shuffle=False)

    return matrix


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression, csq_fields: str
) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns and array of VEP VCF CSQ strings (1 per csq in the struct).

    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.

    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    :param vep_expr: The input VEP Struct
    :param csq_fields: The | delimited list of fields to include in the CSQ (in that order)
    :return: The corresponding CSQ strings
    """
    _csq_fields = [f.lower() for f in csq_fields.split("|")]

    def get_csq_from_struct(
        element: hl.expr.StructExpression, feature_type: str
    ) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                "allele": element.variant_allele,
                "consequence": hl.delimit(element.consequence_terms, delimiter="&"),
                "feature_type": feature_type,
                "feature": (
                    element.transcript_id
                    if "transcript_id" in element
                    else element.regulatory_feature_id
                    if "regulatory_feature_id" in element
                    else element.motif_feature_id
                    if "motif_feature_id" in element
                    else ""
                ),
                "variant_class": vep_expr.variant_class,
            }
        )

        # Only interested in transcripts for now
        fields.update(
            {
                "canonical": hl.cond(element.canonical == 1, "YES", ""),
                "ensp": element.protein_id,
                "gene": element.gene_id,
                "symbol": element.gene_symbol,
                "symbol_source": element.gene_symbol_source,
                "cdna_position": hl.str(element.cdna_start)
                + hl.cond(
                    element.cdna_start == element.cdna_end,
                    "",
                    "-" + hl.str(element.cdna_end),
                ),
                "cds_position": hl.str(element.cds_start)
                + hl.cond(
                    element.cds_start == element.cds_end,
                    "",
                    "-" + hl.str(element.cds_end),
                ),
                "protein_position": hl.str(element.protein_start)
                + hl.cond(
                    element.protein_start == element.protein_end,
                    "",
                    "-" + hl.str(element.protein_end),
                ),
                "sift": element.sift_prediction
                + "("
                + hl.format("%.3f", element.sift_score)
                + ")",
                "polyphen": element.polyphen_prediction
                + "("
                + hl.format("%.3f", element.polyphen_score)
                + ")",
            }
        )

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, "")), "") for f in _csq_fields], "|"
        )

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(
            vep_expr["transcript_consequences"].map(
                lambda x: get_csq_from_struct(x, feature_type="Transcript")
            ),
            hl.empty_array(hl.tstr),
        )
    )

    # prior filtering on consequence will make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def extract_broad_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)

    replace with placeholder if empty
    placeholder values should be least consequential
    e.g. most tools score 0, but for Sift 1 is least important

    :param matrix:
    :return: input matrix with annotations pulled into INFO
    """
    # filter the matrix table on per-consequence basis
    return matrix.annotate_rows(
        info=matrix.info.annotate(
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


def filter_to_classified(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filters to only rows which have an associated class
    :param matrix:
    :return: input matrix, minus rows without Classes applied
    """
    return matrix.filter_rows(
        (matrix.info.Class1 == 1)
        | (matrix.info.Class2 == 1)
        | (matrix.info.Class3 == 1)
        | (matrix.info.Class4 == 1)
    )


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
    # temp - dump this with and without filtering
    logging.info('Filtering Variant rows')
    matrix = filter_mt_rows(
        matrix=matrix, config=config_dict, green_genes=green_gene_set_expression
    )

    # pull annotations into info (not vep)
    logging.info('Pulling VEP annotations into INFO field')
    matrix = extract_broad_annotations(matrix)

    # add Classes to the MT
    logging.info('Applying classes to variant consequences')
    matrix = annotate_class_1(matrix, config_dict)
    matrix = annotate_class_3(matrix, config_dict)
    matrix = annotate_class_2(matrix, config_dict)
    matrix = annotate_class_4(matrix, config_dict)

    # possibly add a background class here for interesting, but only
    # good enough to be a second hit. C4 is this for now

    # filter to class-annotated only prior to export
    logging.info('Filter variants to leave only classified')
    matrix = filter_to_classified(matrix)

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            CSQ=vep_struct_to_csq(matrix.vep, csq_fields=config_dict.get('csq_string'))
        )
    )

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
