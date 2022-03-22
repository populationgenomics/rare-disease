"""
Read, filter, annotate, classify, and write Genetic data
- read VCF into MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- annotate
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- consequence filter
- remove all rows with no interesting consequences
- extract vep data into CSQ string(s)
- annotate with classes 1, 2, 3, and 4
- remove all un-classified variants
- write as VCF

This doesn't include applying inheritance pattern filters
Classifications here are treated as unconfirmed
"""

from typing import Any, Dict, Optional, Tuple
from itertools import permutations
import json
import logging
import sys

import click
import hail as hl

from cloudpathlib import AnyPath


# set some Hail constants
MISSING_STRING = hl.str('missing')
MISSING_INT = hl.int32(0)
ONE_INT = hl.int32(1)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_FLOAT_HI = hl.float64(1.0)

CONFLICTING = hl.str('conflicting')
LOFTEE_HC = hl.str('HC')
PATHOGENIC = hl.str('pathogenic')


def annotate_class_1(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the Class1 annotation (1 or 0) as appropriate
    semi-rare in Gnomad
    at least one Clinvar star
    either Pathogenic or Likely_pathogenic in Clinvar

    Didn't handle 'Pathogenic/Likely_pathogenic'
    Changing to 'contains pathogenic and not conflicting'
    :param matrix:
    :return: same Matrix, with additional field per variant
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class1=hl.if_else(
                (matrix.info.clinvar_stars > 0)
                & (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                & ~(matrix.info.clinvar_sig.lower().contains(CONFLICTING)),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_class_2(
    matrix: hl.MatrixTable, config: Dict[str, Any], new_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    - Gene is new in PanelApp
    - Rare in Gnomad, and
    - Clinvar, or
    - Critical protein consequence on at least one transcript
    - High in silico consequence

    New update! this is now restricted to the NEW genes only
    This means that these are now confident Class2, and we
    only have a MOI test remaining

    :param matrix:
    :param config:
    :param new_genes: the new genes in this panelapp content
    :return: same Matrix, with additional field per variant
    """

    critical_consequences = hl.set(config.get('critical_csq'))

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class2=hl.if_else(
                (new_genes.contains(matrix.geneIds))
                & (
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
                    | (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                    | (
                        (matrix.info.cadd > config['in_silico']['cadd'])
                        | (matrix.info.revel > config['in_silico']['revel'])
                    )
                ),
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
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar
    :param matrix:
    :param config:
    :return:
    """

    critical_consequences = hl.set(config.get('critical_csq'))

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
                & (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: (x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))
                        )
                    )
                    | (matrix.info.clinvar_sig.lower().contains(PATHOGENIC))
                ),
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
    - polyphen and sift are evaluated per-consequence
    :param matrix:
    :param config:
    :return:
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            Class4=hl.if_else(
                (
                    (matrix.info.cadd > config.get('cadd'))
                    & (matrix.info.revel > config.get('revel'))
                )
                | (
                    (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.sift_score, MISSING_FLOAT_HI)
                            < config.get('sift')
                        )
                    )
                    & (
                        matrix.vep.transcript_consequences.any(
                            lambda x: hl.or_else(x.polyphen_score, MISSING_FLOAT_LO)
                            > config.get('polyphen')
                        )
                    )
                    & (
                        (matrix.info.mutationtaster.contains('D'))
                        | (matrix.info.mutationtaster == 'missing')
                    )
                    & (matrix.info.gerp_rs >= config.get('gerp'))
                    & (matrix.info.eigen_phred > config.get('eigen'))
                ),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def annotate_class_4_only(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies a flagto all variants with only Class4 flags applied
    this becomes relevant when we are looking at variants eligible for
    compound het analysis

    :param matrix:
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            class_4_only=hl.if_else(
                (matrix.info.Class1 == 0)
                & (matrix.info.Class2 == 0)
                & (matrix.info.Class3 == 0)
                & (matrix.info.Class4 == 1),
                ONE_INT,
                MISSING_INT,
            )
        )
    )


def extract_comp_het_details(matrix: hl.MatrixTable, write_path: str):
    """
    takes the matrix table, and finds compound-hets per sample
    based on the gene name only

    :param matrix:
    :param write_path: path to the output VCF, transform to json
    """

    # set a new group of values as the key, so that we can collect on them easily
    ch_matrix = matrix.key_rows_by(
        matrix.locus, matrix.alleles, matrix.info.class_4_only
    )
    ch_matrix = ch_matrix.annotate_cols(
        hets=hl.agg.group_by(
            ch_matrix.info.gene,
            hl.agg.filter(ch_matrix.GT.is_het(), hl.agg.collect(ch_matrix.row_key)),
        )
    )

    # extract those possible compound het pairs out as a non-Hail structure
    hets_as_python = ch_matrix.select_cols('hets').col.collect()

    # find all variant pair permutations which aren't both class 4
    combinations = {
        row.s: {
            key: [
                (var1, var2)
                for var1, var2 in permutations(value, 2)
                if not all(
                    map(lambda x: x == 0, [var1['class_4_only'], var2['class_4_only']])
                )
            ]
            for key, value in dict(row.hets).items()
        }
        for row in hets_as_python
    }

    # swap file extensions, to create a JSON output
    output_path = f'{write_path.split(".", maxsplit=1)[0]}.json'

    # write this content to an output path
    with open(AnyPath(output_path), encoding='utf-8') as handle:
        json.dump(combinations, handle)


def filter_matrix_by_ac(
    matrix_data: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """

    :param matrix_data:
    :param config:
    :return:
    """

    # count the samples in the VCF, and use to decide whether to implement
    # 'common within this joint call' as a filter
    # if we reach the sample threshold, filter on AC
    if matrix_data.count_cols() >= config['min_samples_to_ac_filter']:
        matrix_data = matrix_data.filter_rows(
            matrix_data.info.AC <= matrix_data.info.AN // config['ac_filter_percentage']
        )
    return matrix_data


def filter_matrix_by_variant_attributes(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with normalised, high quality variants
    :param matrix_data:
    :return:
    """

    # hard filter for quality; assuming data is well normalised in pipeline
    matrix_data = matrix_data.filter_rows(
        (
            matrix_data.filters.length() == 0
        )  # clarify with FILTER='PASS' & GATK SNP tranches
        & (hl.len(matrix_data.alleles) == 2)
        & (matrix_data.alleles[1] != '*')
    )
    return matrix_data


def annotate_using_vep(matrix_data: hl.MatrixTable) -> hl.MatrixTable:
    """
    runs VEP annotation on the MT
    :param matrix_data:
    :return:
    """

    # throw in a repartition here (annotate even chunks in parallel)
    matrix_data = matrix_data.repartition(150, shuffle=False)

    # VEP 105 annotations
    return hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')


def filter_rows_for_rare(
    matrix: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad & exac
    :param matrix:
    :param config:
    :return:
    """
    # exac and gnomad must be below threshold or missing
    # if missing they were previously replaced with 0.0
    # could also extend this filter to include max gnomad Homs
    return matrix.filter_rows(
        (matrix.info.exac_af < config['af_semi_rare'])
        & (matrix.info.gnomad_af < config['af_semi_rare'])
    )


def filter_benign(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter out benign variants, where clinvar is confident
    :param matrix:
    :return:
    """
    # remove all clinvar benign, decent level of support
    benign = hl.str('benign')
    return matrix.filter_rows(
        (matrix.info.clinvar_sig.lower().contains(benign))
        & (matrix.info.clinvar_stars > 0),
        keep=False,
    )


def filter_to_green_genes_and_split(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    reduces geneIds set to green only, then splits
    :param matrix:
    :param green_genes:
    :return:
    """

    # remove any rows with no genic annotation at all
    # do this here as set intersections with Missing will fail
    matrix = matrix.filter_rows(hl.is_missing(matrix.geneIds), keep=False)

    # replace the default list of green IDs with a reduced set
    matrix = matrix.annotate_rows(geneIds=green_genes.intersection(matrix.geneIds))

    # split to form a separate row for each green gene
    # this transforms the 'geneIds' field from a set to a string
    return matrix.explode_rows(matrix.geneIds)


def filter_by_consequence(
    matrix: hl.MatrixTable, config: Dict[str, Any]
) -> hl.MatrixTable:
    """
    - reduce the per-row transcript consequences to those specific to the geneIds
    - reduce the rows to ones where there are remaining tx consequences

    :param matrix:
    :param config: dictionary content relating to hail
    :return: reduced matrix
    """

    # identify consequences to discard from the config
    useless_csq = hl.set(config['useless_csq'])

    # reduce consequences to overlap with per-variant green geneIDs (pre-filtered)
    # added another condition to state that the tx biotype needs to be protein_coding,
    # unless the row also has an attached MANE transcript
    # consider an extra allowance for strong Appris transcripts
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: (matrix.geneIds == x.gene_id)
                & (hl.len(hl.set(x.consequence_terms).difference(useless_csq)) > 0)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    # filter out all rows with no remaining consequences
    return matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) > 0)


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression, csq_fields: str
) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns an array of VEP VCF CSQ strings
    (1 per csq in the struct).

    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.

    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    :param vep_expr: The input VEP Struct
    :param csq_fields: | delimited fields to include in the CSQ (in that order)
    :return: The corresponding CSQ strings
    """
    _csq_fields = [f.lower() for f in csq_fields.split('|')]

    def get_csq_from_struct(
        element: hl.expr.StructExpression, feature_type: str
    ) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'allele': element.variant_allele,
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature_type': feature_type,
                'feature': (
                    element.transcript_id
                    if 'transcript_id' in element
                    else element.regulatory_feature_id
                    if 'regulatory_feature_id' in element
                    else element.motif_feature_id
                    if 'motif_feature_id' in element
                    else ''
                ),
                'variant_class': vep_expr.variant_class,
                'canonical': hl.if_else(element.canonical == 1, 'YES', ''),
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'symbol_source': element.gene_symbol_source,
                'cdna_position': hl.str(element.cdna_start)
                + hl.if_else(
                    element.cdna_start == element.cdna_end,
                    '',
                    '-' + hl.str(element.cdna_end),
                ),
                'cds_position': hl.str(element.cds_start)
                + hl.if_else(
                    element.cds_start == element.cds_end,
                    '',
                    '-' + hl.str(element.cds_end),
                ),
                'protein_position': hl.str(element.protein_start)
                + hl.if_else(
                    element.protein_start == element.protein_end,
                    '',
                    '-' + hl.str(element.protein_end),
                ),
                'sift': element.sift_prediction
                + '('
                + hl.format('%.3f', element.sift_score)
                + ')',
                'polyphen': element.polyphen_prediction
                + '('
                + hl.format('%.3f', element.polyphen_score)
                + ')',
                'mane_select': element.mane_select,
            }
        )

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, '')), '') for f in _csq_fields], '|'
        )

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(
            vep_expr['transcript_consequences'].map(
                lambda x: get_csq_from_struct(x, feature_type='Transcript')
            ),
            hl.empty_array(hl.tstr),
        )
    )

    # prior filtering on consequence will make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def extract_annotations(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)

    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is least important

    :param matrix:
    :return: input matrix with annotations pulled into INFO
    """

    return matrix.annotate_rows(
        info=matrix.info.annotate(
            exac_af=hl.or_else(matrix.exac.AF, MISSING_FLOAT_LO),
            exac_ac_het=hl.or_else(matrix.exac.AC_Het, MISSING_INT),
            exac_ac_hom=hl.or_else(matrix.exac.AC_Hom, MISSING_INT),
            exac_ac_hemi=hl.or_else(matrix.exac.AC_Hemi, MISSING_INT),
            gnomad_ex_cov=hl.or_else(matrix.gnomad_exome_coverage, MISSING_FLOAT_LO),
            gnomad_ex_af=hl.or_else(matrix.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_ex_an=hl.or_else(matrix.gnomad_exomes.AN, MISSING_INT),
            gnomad_ex_ac=hl.or_else(matrix.gnomad_exomes.AC, MISSING_INT),
            gnomad_ex_hom=hl.or_else(matrix.gnomad_exomes.Hom, MISSING_INT),
            gnomad_cov=hl.or_else(matrix.gnomad_genome_coverage, MISSING_FLOAT_LO),
            gnomad_af=hl.or_else(matrix.gnomad_genomes.AF, MISSING_FLOAT_LO),
            gnomad_an=hl.or_else(matrix.gnomad_genomes.AN, MISSING_INT),
            gnomad_ac=hl.or_else(matrix.gnomad_genomes.AC, MISSING_INT),
            gnomad_hom=hl.or_else(matrix.gnomad_genomes.Hom, MISSING_INT),
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
            # these next 3 are per-transcript, with ';' to delimit
            # pulling these annotations into INFO with ';' to separate
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
    Filter to rows tagged with a class
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
    write the remaining MatrixTable content to file as a VCF
    :param matrix:
    :param output_path: where to write
    """
    hl.export_vcf(
        matrix,
        output_path,
        tabix=True,
    )


def green_and_new_from_panelapp(
    panel_data: Dict[str, Dict[str, str]]
) -> Tuple[hl.SetExpression, hl.SetExpression]:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes
    Also identify the subset of those genes which relate to NEW in panel
    :param panel_data:
    :return: two set expressions, green genes and new genes
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_data.keys()) - {'panel_metadata'}
    logging.info('Extracted %d green genes', len(green_genes))
    green_gene_set_expression = hl.literal(green_genes)

    new_genes = {gene for gene in green_genes if panel_data[gene].get('new')}
    logging.info('Extracted %d NEW genes', len(new_genes))
    new_gene_set_expression = hl.literal(new_genes)

    return green_gene_set_expression, new_gene_set_expression


@click.command()
@click.option('--mt', 'mt_path', help='path to the matrix table to ingest')
@click.option('--pap', 'panelapp_path', help='bucket path containing panelapp JSON')
@click.option('--config', 'config_path', help='path to a config dict')
@click.option('--output', 'out_vcf', help='VCF path to export results')
@click.option(
    '--mt_out',
    help='path to export annotated MT to',
    required=False,
    default=None,
)
def main(
    mt_path: str,
    panelapp_path: str,
    config_path: str,
    out_vcf: str,
    mt_out: Optional[str] = None,
):
    """
    Read the MT from disk
    Do filtering and class annotation
    Export as a VCF

    :param mt_path: path to the MT directory
    :param panelapp_path: path to the panelapp data dump
    :param config_path: path to the config json
    :param out_vcf: path to write the VCF out to
    :param mt_out:
    """

    # get the run configuration JSON
    logging.info('Reading config dict from "%s"', config_path)
    with open(AnyPath(config_path), encoding='utf-8') as handle:
        config_dict = json.load(handle)

    # find the config area specific to hail operations
    hail_config = config_dict.get('hail')

    # read the parsed panelapp data
    logging.info('Reading PanelApp data from "%s"', panelapp_path)
    with open(AnyPath(panelapp_path), encoding='utf-8') as handle:
        panelapp = json.load(handle)

    # pull green and new genes from the panelapp data
    green_expression, new_expression = green_and_new_from_panelapp(panelapp)

    logging.info(
        'Starting Hail with reference genome "%s"', hail_config.get('ref_genome')
    )
    # initiate Hail with the specified reference
    hl.init(default_reference=hail_config.get('ref_genome'), quiet=True)

    # if we already generated the annotated output, load instead
    if AnyPath(mt_out.rstrip('/') + '/').exists():
        logging.info('Loading annotated MT from "%s"', mt_out)
        matrix = hl.read_matrix_table(mt_out)
        logging.debug('Loaded annotated MT, size: %d', matrix.count_rows())

    else:
        logging.info('Loading MT from "%s"', mt_path)
        matrix = hl.read_matrix_table(mt_path)
        logging.debug('Loaded new MT, size: %d', matrix.count_rows())

        # hard filter entries in the MT prior to annotation
        logging.info('Hard filtering variants')
        matrix = filter_matrix_by_ac(matrix_data=matrix, config=hail_config)
        matrix = filter_matrix_by_variant_attributes(matrix_data=matrix)

        # re-annotate using VEP
        logging.info('Annotating variants')
        matrix = annotate_using_vep(matrix_data=matrix)

        # if a path is provided, dump the MT
        if mt_out is not None:
            matrix.write(mt_out, overwrite=True)

    # pull annotations into info and update if missing
    logging.info('Pulling VEP annotations into INFO field')
    matrix = extract_annotations(matrix)

    # filter on row annotations
    logging.info('Filtering Variant rows')
    matrix = filter_rows_for_rare(matrix=matrix, config=hail_config)
    matrix = filter_benign(matrix=matrix)
    matrix = filter_to_green_genes_and_split(
        matrix=matrix, green_genes=green_expression
    )
    matrix = filter_by_consequence(matrix=matrix, config=hail_config)

    logging.info('Repartition to 50 fragments following Gene ID filter')
    matrix = matrix.repartition(50, shuffle=False)

    # add Classes to the MT
    logging.info('Applying classes to variant consequences')
    matrix = annotate_class_1(matrix)
    matrix = annotate_class_2(matrix, hail_config, new_expression)
    matrix = annotate_class_3(matrix, hail_config)
    matrix = annotate_class_4(matrix, hail_config['in_silico'])

    # filter to class-annotated only prior to export
    logging.info('Filter variants to leave only classified')
    matrix = filter_to_classified(matrix)

    # add an additional annotation, if the variant is Class4 only
    matrix = annotate_class_4_only(matrix)

    # parse out the compound het details
    extract_comp_het_details(matrix=matrix, write_path=out_vcf)

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    matrix = matrix.annotate_rows(
        info=matrix.info.annotate(
            CSQ=vep_struct_to_csq(
                matrix.vep, csq_fields=config_dict['variant_object'].get('csq_string')
            ),
            gene_id=matrix.geneIds,
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
