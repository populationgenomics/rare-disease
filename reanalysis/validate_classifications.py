"""
runs between classification and publishing results
takes 2 VCFs: classes and compound hets
reads in all compound het pairs
reads in all panelapp details
for each variant in each participant, check MOI
"""
import logging
from typing import Any, Dict, List, Union
from itertools import chain
import json
import os
import click
from cyvcf2 import VCFReader
from google.cloud import storage

from reanalysis.moi_tests import MOIRunner
from reanalysis.results_builder import HTMLBuilder
from reanalysis.utils import (
    AnalysisVariant,
    VARIANT_STRING_TEMPLATE,
    COMP_HET_VALUES,
    get_simple_moi,
    parse_ped_simple,
    PedPerson,
    ReportedVariant,
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
            paired_var_string = VARIANT_STRING_TEMPLATE.format(
                row_dict.get('chrom').replace('chr', ''),
                row_dict.get('pos'),
                row_dict.get('ref'),
                row_dict.get('alt'),
            )

            # don't allow a 'pair' with the current variant
            if (
                int(row_dict['pos']) == variant.POS
                and row_dict['ref'] == variant.REF
                and row_dict['alt'] == variant.ALT[0]
            ):
                # logging.error(
                #     'Found a variant forming a comp-het with itself: %s',
                #     paired_var_string,
                # )
                continue

            # noinspection PyTypeChecker
            comp_het_lookup.setdefault(row_dict.get('sample'), {})[
                f'{paired_var_string}'
            ] = AnalysisVariant(variant, samples)

    logging.info('%d samples contain a compound het', len(comp_het_lookup))
    # for sample in comp_het_lookup.values():
    #     for variant in sample.values():
    #         print(repr(variant.var), variant.get_class_ints())

    # noinspection PyTypeChecker
    return comp_het_lookup


def set_up_inheritance_filters(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
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
                pedigree=pedigree, target_moi=gene_moi, ad_threshold=ad_threshold
            )

    return moi_dictionary


def validate_class_2(
    panelapp_data: Dict[str, Dict[str, Union[str, bool]]],
    variant: AnalysisVariant,
    moi_lookup: Dict[str, MOIRunner],
    comp_het_lookup: Dict[str, Dict[str, AnalysisVariant]],
) -> bool:
    """
    Test for class 2 variant
    Is it new in PanelApp?
    Otherwise, does it pass more MOI tests than it did previously?
    :param panelapp_data:
    :param variant:
    :param moi_lookup:
    :param comp_het_lookup:
    :return: boolean, whether to bother analysing
    """

    # start with the assumption we will discard
    retain = False

    gene = variant.var.INFO.get('gene_id')

    # get relevant panelapp contents
    panel_gene_data = panelapp_data.get(gene)

    # if variant is C2, we need the time sensitive check
    # get the panelapp differences
    # might need some more digging here to clarify logic **
    if panel_gene_data.get('new'):
        retain = True

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
) -> List[ReportedVariant]:
    """

    :param classified_variant_source:
    :param comp_het_lookup:
    :param moi_lookup:
    :param panelapp_data:
    :return:
    """

    results = []

    variant_source = VCFReader(classified_variant_source)
    vcf_samples = variant_source.samples

    # crawl through all the results
    for variant in variant_source:

        # cast as an analysis variant
        analysis_variant = AnalysisVariant(variant, samples=vcf_samples)

        gene = variant.INFO.get('gene_id')

        # #  allow for multiple genes here?
        # # e.g. iterate over all possible genes
        # # store class 2 status, check for _this gene_, then reset at the start of the loop
        # c2_status = analysis_variant.class_2
        # for gene_id in gene:
        #     analysis_variant.class_2 = c2_status
        #     # then do C2 test, and follow with MOI test

        # one variant appears to be retained here, in a red gene
        # possibly overlapping with a Green gene?
        if gene not in panelapp_data:
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
                panelapp_data=panelapp_data,
                moi_lookup=moi_lookup,
                variant=analysis_variant,
                comp_het_lookup=comp_het_lookup,
            ):
                analysis_variant.class_2 = False

        # we never use a C4-only variant as a principal variant
        # and we don't consider a variant with no assigned classes
        # if not analysis_variant.is_classified():
        if analysis_variant.class_4_only() or not analysis_variant.is_classified():
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
    There was a possibility that a single variant can be
    classified in different ways. This method cleans those
    down to unique

    Join all possible classes for the condensed variants
    :param result_list:
    :return:
    """

    clean_results: Dict[str, Dict[str, ReportedVariant]] = {}

    for each_instance in result_list:
        support_id = (
            repr(each_instance.support_var.var)
            if each_instance.support_var is not None
            else 'Unsupported'
        )
        var_uid = (
            f'{repr(each_instance.var_data.var)}__'
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
    results = apply_moi_to_variants(
        classified_variant_source=classified_vcf,
        comp_het_lookup=comp_het_digest,
        moi_lookup=moi_lookup,
        panelapp_data=panelapp_data,
    )

    # remove duplicates of the same variant
    cleaned_results = clean_initial_results(results)

    # use the config file to select the relevant CPG to Seqr ID JSON file
    seqr_data = read_json_dict_from_path(config_dict.get('seqr_lookup'))

    # generate some html
    html_maker = HTMLBuilder(
        results_dict=cleaned_results, seqr_lookup=seqr_data, panelapp_data=panelapp_data
    )
    html_tables, class_2_genes = html_maker.create_html_tables()

    class_2_table = html_maker.class_2_table(class_2_genes)

    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write('<head>\n</head>\n<body>\n')
        handle.write('<h3>MOI changes used for Class 2</h3>')
        handle.write(class_2_table)
        handle.write('<br/>')

        handle.write('<h1>Per Sample Results</h1>')
        for sample, table in html_tables.items():
            handle.write(fr'<h3>{sample}</h3>')
            handle.write(table)
        handle.write('\n</body>')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
