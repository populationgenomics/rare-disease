"""
runs between classification and publishing results
takes 2 VCFs: classes and compound hets
reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI
"""
import json
import logging
from typing import Any, Dict, List, Union
from itertools import chain
import click
from cyvcf2 import VCFReader

from reanalysis.moi_tests import MOIRunner
from reanalysis.results_builder import HTMLBuilder
from reanalysis.utils import (
    AnalysisVariant,
    COMP_HET_VALUES,
    CustomEncoder,
    get_simple_moi,
    parse_ped_simple,
    PedPerson,
    read_json_dictionary,
    ReportedVariant,
    VARIANT_STRING_TEMPLATE,
)


def parse_comp_hets(
    vcf_path: str, config: Dict[str, Any]
) -> Dict[str, Dict[str, AnalysisVariant]]:
    """
    iterate over variants, and associate samples with variant pairs
    check source code, but slivar should only pair on same transcript
    {
      sample: {
        'chr-pos-ref-alt': paired_variant,
      }
    }

    :param vcf_path:
    :param config:
    """
    comp_het_lookup = {}

    var_source = VCFReader(vcf_path)

    samples = var_source.samples

    # iterate over all variants
    for variant in var_source:

        # create a variant object instance
        a_variant = AnalysisVariant(variant, samples, config=config)

        # get all comp-het annotations matching this position
        # expect one, allow for surprises
        for var_info in a_variant.info.get('slivar_comphet').split(','):

            # cast info into a dictionary
            row_dict = dict(zip(COMP_HET_VALUES, var_info.split('/')))

            # create a string
            paired_var_string = VARIANT_STRING_TEMPLATE.format(
                row_dict.get('chrom').replace('chr', ''),
                row_dict.get('pos'),
                row_dict.get('ref'),
                row_dict.get('alt'),
            )

            # don't allow a 'pair' with the current variant
            if paired_var_string == a_variant.string:
                continue

            # noinspection PyTypeChecker
            comp_het_lookup.setdefault(row_dict.get('sample'), {})[
                f'{paired_var_string}'
            ] = a_variant

    logging.info('%d samples contain a compound het', len(comp_het_lookup))

    # noinspection PyTypeChecker
    return comp_het_lookup


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    config: Dict[str, Any],
    pedigree: Dict[str, PedPerson],
) -> Dict[str, MOIRunner]:
    """
    parse the panelapp data, and find all MOIs in this dataset
    for each unique MOI, set up a MOI filter instance
    save each one to a dictionary

    {MOI_string: MOI_runner (with a .run() method)}

    The MOI_runner class will use the provided MOI string to
    select which filters will be appropriate

    All logic regarding how MOI is applied, and which MOIs to
    apply to which PanelApp MOI descriptions is partitioned off into
    the MOI class. All we need here is a Run() method, that returns
    either a list of results, or an empty list

    for every variant, we can then do a simple lookup using this
    dictionary to find the correct MOI runner, and run it
    that will return all matching MOIs for the variant

    This dictionary format means we only have to set up each once
    A billion variants, 6 MOI = 6 test instances, each created once
    :param panelapp_data:
    :param config:
    :param pedigree:
    :return:
    """

    moi_dictionary = {}

    # iterate over all genes
    for key, gene_data in panelapp_data.items():

        # skip over the stored metadata
        if '_version' in key:
            continue

        # extract the per-gene MOI, and SIMPLIFY
        gene_moi = get_simple_moi(gene_data.get('moi'))

        # if we haven't seen this MOI before, set up the appropriate filter
        if gene_moi not in moi_dictionary:

            # get a MOIRunner with the relevant filters
            moi_dictionary[gene_moi] = MOIRunner(
                pedigree=pedigree, target_moi=gene_moi, config=config['moi_tests']
            )

    return moi_dictionary


def validate_class_2(
    gene: str,
    panel_gene_data: Dict[str, str],
    variant: AnalysisVariant,
    moi_lookup: Dict[str, MOIRunner],
    comp_het_lookup: Dict[str, Dict[str, AnalysisVariant]],
    new_only: bool = False,
) -> bool:
    """
    Test for class 2 variant
    Is it new in PanelApp?
    Otherwise, does it pass more MOI tests than it did previously?

    ** this latter test is currently not being used
    The comparison will initially be
    "is the variant in a gene which is green in the Mendeliome now,
    but was not in the pre-PanelApp Mendeliome"


    :param gene:
    :param panel_gene_data:
    :param variant:
    :param moi_lookup:
    :param comp_het_lookup:
    :param new_only: if true, only use newly green values in a Class 2 test
    :return: boolean, whether to bother analysing
    """

    # start with the assumption we will discard
    retain = False

    # if variant is C2, we need the time sensitive check
    # get the panelapp differences
    # might need some more digging here to clarify logic **
    if panel_gene_data.get('new'):
        retain = True

    # if we're only interested in fully NEW green genes, skip MOI tests
    elif new_only:
        retain = False

    # check there are more classifying events now than previously
    elif panel_gene_data.get('changed'):

        # the 'simplified' MOI for 'Unknown' is 'search for both mono and biallelic inheritance'
        # if we use the simplified version, we fail to identify older MOI being unpopulated
        # e.g, Unknown -> Monoallelic
        # this situation could actually return more results as "Unknown", given that it's
        # evaluated as Dominant and Recessive
        if panel_gene_data.get('old_moi') == 'Unknown':
            retain = True

        else:
            old_moi = get_simple_moi(panel_gene_data.get('old_moi'))
            new_moi = get_simple_moi(panel_gene_data.get('moi'))

            # get successful tiering modes for the new and old MOIs
            # collected into a set of Strings
            old_moi_passes = set(
                chain.from_iterable(
                    [
                        result.reasons
                        for result in moi_lookup[old_moi].run(
                            principal_var=variant, comp_hets=comp_het_lookup, ensg=gene
                        )
                    ]
                )
            )
            new_moi_passes = set(
                chain.from_iterable(
                    [
                        result.reasons
                        for result in moi_lookup[new_moi].run(
                            principal_var=variant, comp_hets=comp_het_lookup, ensg=gene
                        )
                    ]
                )
            )

            # require at least one using current MOI,
            # and a difference between old and new
            retain = len(new_moi_passes - old_moi_passes) > 0

    return retain


# pylint: disable=too-many-locals
def apply_moi_to_variants(
    classified_variant_source: str,
    comp_het_lookup: Dict[str, Dict[str, AnalysisVariant]],
    moi_lookup: Dict[str, MOIRunner],
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    config: Dict[str, Any],
) -> List[ReportedVariant]:
    """

    :param classified_variant_source:
    :param comp_het_lookup:
    :param moi_lookup:
    :param panelapp_data:
    :param config:
    :return:
    """

    # parse from config... should we use NEW only for Cat. 2?
    # alternative is considering both NEW and changed MOI
    class_2_new_only = config['moi_tests'].get('class_2_new_only', True)

    results = []

    variant_source = VCFReader(classified_variant_source)
    vcf_samples = variant_source.samples

    # crawl through all the results
    for variant in variant_source:

        # cast as an analysis variant
        analysis_variant = AnalysisVariant(variant, samples=vcf_samples, config=config)

        # each variant is uniformly associated with a single gene
        # each variant also has at least one tx_csq, so take the first
        # CHANGE - take geneIds, copy into INFO, split on that field instead
        gene = analysis_variant.info['transcript_consequences'][0].get('gene')

        # extract the panel data specific to this gene
        panel_gene_data = panelapp_data.get(gene)

        # one variant appears to be retained here, in a red gene
        # possibly overlapping with a Green gene?
        if panel_gene_data is None:
            logging.error("How did this gene creep in? %s", gene)
            continue

        # if variant is C2, we need the time sensitive check
        # get the panelapp differences
        if analysis_variant.class_2:

            # if there was only one panel version, definitely skip
            # if we don't have a _new_ MOI to use, or a new gene, skip
            if panelapp_data["panel_metadata"].get(
                "previous_version"
            ) is None or not validate_class_2(
                gene=gene,
                panel_gene_data=panel_gene_data,
                moi_lookup=moi_lookup,
                variant=analysis_variant,
                comp_het_lookup=comp_het_lookup,
                new_only=class_2_new_only,
            ):
                analysis_variant.class_2 = False

        # we never use a C4-only variant as a principal variant
        # and we don't consider a variant with no assigned classes
        if analysis_variant.class_4_only or not analysis_variant.is_classified:
            continue

        results.extend(
            moi_lookup[get_simple_moi(panelapp_data[gene].get('moi'))].run(
                principal_var=analysis_variant, comp_hets=comp_het_lookup, ensg=gene
            )
        )

    return results


def clean_initial_results(
    result_list: List[ReportedVariant],
) -> Dict[str, Dict[str, ReportedVariant]]:
    """
    Possibility 1 variant can be classified multiple ways
    This cleans those to unique for final report
    Join all possible classes for the condensed variants
    :param result_list:
    """

    clean_results: Dict[str, Dict[str, ReportedVariant]] = {}

    for each_instance in result_list:
        support_id = (
            'Unsupported'
            if not each_instance.supported
            else each_instance.support_var.string
        )
        var_uid = (
            f'{each_instance.var_data.string}__'
            f'{each_instance.gene}__'
            f'{support_id}'
        )

        # create a section for this sample if it doesn't exist
        if each_instance.sample not in clean_results:
            clean_results[each_instance.sample] = {}

        # get an existing object, or use the current one
        variant_object = clean_results.setdefault(each_instance.sample, {}).setdefault(
            var_uid, each_instance
        )

        # combine any possible reasons, and add
        clean_results[each_instance.sample][
            var_uid
        ].reasons = variant_object.reasons.union(each_instance.reasons)
    return clean_results


@click.command()
@click.option('--conf', 'config_path', help='Path to a config JSON file')
@click.option(
    '--class_vcf',
    help='VCF from Hail with variant categories and annotations',
)
@click.option('--comp_het', help='VCF limited to compound-hets')
@click.option('--ped', 'pedigree', help='Pedigree file')
@click.option('--pap', 'panelapp', help='PanelApp JSON file')
@click.option('--out_path', help='Where to write the output to')
@click.option('--out_json', help='Write the analysis results in JSON form')
def main(
    config_path: str,
    class_vcf: str,
    comp_het: str,
    pedigree: str,
    panelapp: str,
    out_path: str,
    out_json: str,
):  # pylint: disable=too-many-arguments
    """
    All VCFs in use at this point will be small
    These have been pre-filtered to retain only a small number of classified variants
    holding all the variants in memory should not be a challenge, no matter how large
    the cohort; if the variant number is large, the classes should be refined
    We expect approximately linear scaling with participants in the joint call

    Might be able to use a single output path, just altering the extension
    Depends on how this is handled by Hail, as the object paths are Resource File paths

    :param config_path:
    :param class_vcf:
    :param comp_het:
    :param pedigree:
    :param panelapp:
    :param out_path:
    :param out_json:
    """

    # parse the pedigree from the file
    pedigree_digest = parse_ped_simple(pedigree)

    # parse panelapp data from dict
    panelapp_data = read_json_dictionary(panelapp)

    # get the runtime configuration
    config_dict = read_json_dictionary(config_path)

    # find all the Compound Hets from CH VCF
    comp_het_digest = parse_comp_hets(comp_het, config=config_dict)

    # set up the inheritance checks
    moi_lookup = set_up_inheritance_filters(
        panelapp_data=panelapp_data, pedigree=pedigree_digest, config=config_dict
    )

    # find classification events
    results = apply_moi_to_variants(
        classified_variant_source=class_vcf,
        comp_het_lookup=comp_het_digest,
        moi_lookup=moi_lookup,
        panelapp_data=panelapp_data,
        config=config_dict,
    )

    # remove duplicate variants
    cleaned_results = clean_initial_results(results)

    # dump the JSON-friendly results to a file
    with open(out_json, 'w', encoding='utf-8') as handle:
        json.dump(cleaned_results, handle, cls=CustomEncoder, indent=4)

    # generate some html
    html_maker = HTMLBuilder(
        results_dict=cleaned_results,
        panelapp_data=panelapp_data,
        pedigree=pedigree_digest,
        config=config_dict,
    )
    html_maker.write_html(output_path=out_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
