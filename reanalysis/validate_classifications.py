"""
runs between classification and publishing results
takes 2 VCFs: classes and compound hets
reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI
"""
import logging
from typing import Any, Dict
import json
import os
import click
from cyvcf2 import VCFReader
from reanalysis.utils import (
    AnalysisVariant,
    class_2_variant,
    c4_only,
    COMP_HET_TEMPLATE,
    COMP_HET_VALUES,
    get_class_1_2_3_4,
    parse_ped_simple,
    PedPerson,
    string_format_variant,
)
from reanalysis.moi_tests import MOIRunner


def read_json_dictionary(json_path: str) -> Any:
    """
    point at a json file, return the contents
    :param json_path:
    :return: whatever is in the file
    """
    assert os.path.exists(json_path)
    with open(json_path, 'r', encoding='utf-8') as handle:
        json_content = json.load(handle)
    return json_content


def parse_comp_hets(vcf_path: str) -> Dict[str, Dict[str, AnalysisVariant]]:
    """
    iterate over variants, and associate samples with variant pairs
    check source code, but slivar should only pair on same transcript
    {
      sample: {
        'chr-pos-ref-alt-transcript': paired_variant,
      }
    }

    :param vcf_path:
    :return:
    """
    comp_het_lookup = {}

    var_source = VCFReader(vcf_path)

    samples = var_source.samples

    # iterate over all variants
    for variant in var_source:

        # get all comp-het annotations matching this position
        # expect one, allow for surprises
        for var_info in variant.INFO.get('slivar_comphet').split(','):

            # cast info into a dictionary
            row_dict = dict(zip(COMP_HET_VALUES, var_info.split('/')))

            # create a string
            paired_var_string = COMP_HET_TEMPLATE.format(
                row_dict.get('chrom').replace('chr', ''),
                row_dict.get('pos'),
                row_dict.get('ref'),
                row_dict.get('alt'),
                variant.INFO.get('transcript_id'),
            )

            # noinspection PyTypeChecker
            comp_het_lookup.setdefault(row_dict.get('sample'), {})[
                f'{paired_var_string}'
            ] = AnalysisVariant(variant, samples)

    # noinspection PyTypeChecker
    return comp_het_lookup


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Dict[str, str]]],
    config: Dict[str, Any],
    pedigree: Dict[str, PedPerson],
) -> Dict[str, MOIRunner]:
    """

    :param panelapp_data:
    :param config:
    :param pedigree:
    :return:
    """

    moi_dictionary = {}

    # get the stringent threshold to use for dominant MOI
    ad_threshold = config.get('gnomad_dominant')

    # iterate over all genes
    for gene_data in panelapp_data['panel_data'].values():

        # extract the per-gene MOI
        gene_moi = gene_data.get('moi')

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:

            # use the MOIRunner to apply relevant filters
            moi_dictionary[gene_moi] = MOIRunner(
                pedigree=pedigree, target_moi=gene_moi, ad_threshold=ad_threshold
            )

    return moi_dictionary


# pylint: disable=too-many-locals
def apply_moi_to_variants(
    classified_variant_source: str,
    comp_het_lookup: Dict[str, Dict[str, AnalysisVariant]],
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Dict[str, str]]],
):
    """

    :param classified_variant_source:
    :param comp_het_lookup:
    :param moi_lookup:
    :param panelapp_data:
    :return:
    """

    results = []
    # sample_results = []

    variant_source = VCFReader(classified_variant_source)
    vcf_samples = variant_source.samples
    for variant in variant_source:

        # we never use a C4-only variant as a principal variant
        if c4_only(variant):
            continue

        gene = variant.INFO.get('gene_id')

        # cast as an analysis variant
        analysis_variant = AnalysisVariant(variant, samples=vcf_samples)

        # if variant is C2, we need the time sensitive check
        # get the panelapp differences
        # might need some more digging here to clarify logic **
        if class_2_variant(variant):
            if gene in panelapp_data['changes'].get('new'):
                pass
            elif gene in panelapp_data['changes'].get('changed'):
                old_moi, new_moi = panelapp_data['changes'].get('changed').get(gene)

                # consistent - either both passed or failed
                if len(
                    moi_lookup[old_moi].run(
                        principal_var=analysis_variant, comp_hets=comp_het_lookup
                    )
                ) == len(
                    moi_lookup[new_moi].run(
                        principal_var=analysis_variant, comp_hets=comp_het_lookup
                    )
                ):
                    continue
            else:
                continue

        if gene not in panelapp_data['panel_data']:
            logging.error("How did this gene creep in? %s", gene)
            continue
        moi = panelapp_data['panel_data'][gene].get('moi')
        for sample, reason, _variants in moi_lookup[moi].run(
            principal_var=analysis_variant, comp_hets=comp_het_lookup
        ):
            print(
                gene,
                sample,
                reason,
                get_class_1_2_3_4(variant),
                string_format_variant(variant),
            )
    return results


@click.command()
@click.option('--conf', 'config_path', help='')
@click.option('--class_vcf', 'classified_vcf', help='')
@click.option('--comp_het', 'compound_het', help='VCF limited to compound-hets')
@click.option('--ped', 'pedigree', help='Pedigree file')
@click.option('--pap', 'panelapp', help='PanelApp JSON file')
@click.option('--out_path', 'out_path', help='Where to write the output to')
def main(
    config_path: str,
    classified_vcf: str,
    compound_het: str,
    pedigree: str,
    panelapp: str,
    out_path: str,
):  # pylint: disable=too-many-arguments
    """
    All VCFs in use at this point will be small
    These have been pre-filtered to retain only a small number of classified variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call
    :param config_path:
    :param classified_vcf:
    :param compound_het:
    :param pedigree:
    :param panelapp:
    :param out_path:
    :return:
    """

    # parse the pedigree from the file
    pedigree_digest = parse_ped_simple(pedigree)

    # find all the Compound Hets from CH VCF
    comp_het_digest = parse_comp_hets(compound_het)

    # parse panelapp data from dict
    panelapp_data = read_json_dictionary(panelapp)

    # get the runtime configuration
    config_dict = read_json_dictionary(config_path)

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data, pedigree=pedigree_digest, config=config_dict
    )

    # find classification events
    _results = apply_moi_to_variants(
        classified_variant_source=classified_vcf,
        comp_het_lookup=comp_het_digest,
        moi_lookup=moi_lookup,
        panelapp_data=panelapp_data,
    )
    print(out_path)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
