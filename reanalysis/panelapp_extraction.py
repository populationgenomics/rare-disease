"""
PanelApp Parser for Reanalysis

Takes a panel ID
For the latest content, pulls Symbol, ENSG, and MOI
    (MOI is simplified from PanelApp enum)

Optionally user can provide a date in the past
Identify the highest panel version prior to that date
Pull all details from the earlier version
Store all discrepancies between earlier and current

Write all output to a JSON dictionary
"""


import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Union
import json

import requests

import click


# panelapp URL constants
PANELAPP_ROOT = 'https://panelapp.agha.umccr.org/api/v1'
PANEL_ROOT = f'{PANELAPP_ROOT}/panels'
PANEL_CONTENT = f'{PANEL_ROOT}/{{panel_id}}'
ACTIVITIES = f'{PANEL_CONTENT}/activities'


def get_json_response(url: str) -> Union[List[Dict[str, str]], Dict[str, Any]]:
    """
    takes a request URL, checks for healthy response, returns the JSON
    :param url:
    :return:
    """

    response = requests.get(url, headers={'Accept': 'application/json'})
    response.raise_for_status()
    return response.json()


def get_previous_version(panel_id: str, since: datetime) -> str:
    """
    work through the list of panel updates in reverse (earliest -> latest)
    return the first panel version after the threshold date
    If we don't find any, return the latest version

    Note: Django dates include Hours and Minutes, so any Django date from
    the same day as a 'YYYY-MM-DD' date will be 'greater than'
    for the purposes of a value comparison

    :param panel_id: panel ID to use
    :param since: date of the
    :return:
    """

    activity_list = get_json_response(url=ACTIVITIES.format(panel_id=panel_id))
    entry_version = None
    for entry in reversed(activity_list):
        # take the panel version
        entry_version = entry.get('panel_version')
        # uses the django datestamp format
        entry_date = datetime.strptime(entry.get('created'), "%Y-%m-%dT%H:%M:%S.%fz")

        if entry_date >= since:
            return entry_version
    return entry_version


def get_simple_moi(moi: str) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a reduced
    range of cases which can be easily implemented in RD analysis
    :param moi: full PanelApp string
    :return: a simplified representation
    """

    # default to considering both
    panel_app_moi = 'Mono_And_Biallelic'
    if moi is None:
        # exit iteration, return both (all considered)
        return panel_app_moi

    # ideal for match-case, coming to a python 3.10 near you!
    if moi.startswith('BIALLELIC'):
        panel_app_moi = 'Biallelic'
    if moi.startswith("BOTH"):
        panel_app_moi = 'Mono_And_Biallelic'
    if moi.startswith('MONO'):
        panel_app_moi = 'Monoallelic'
    if moi.startswith('X-LINKED'):
        if 'biallelic' in moi:
            panel_app_moi = 'Hemi_Bi_In_Female'
        else:
            panel_app_moi = 'Hemi_Mono_In_Female'

    return panel_app_moi


def get_panel_green(
    reference_genome: Optional[str] = "GRch38",
    panel_id: str = '137',
    version: Optional[str] = None,
) -> Dict[str, Dict[str, Union[datetime, str]]]:
    """
    Takes a panel number, and pulls all details from PanelApp
    :param reference_genome: GRch37 or GRch38
    :param panel_id: defaults to the PanelAppAU Mendeliome
    :param version:
    :return:
    """

    gene_dict = {}

    panel_app_genes_url = PANEL_CONTENT.format(panel_id=panel_id)
    if version is not None:
        panel_app_genes_url += f"?version={version}"

    panel_response = requests.get(panel_app_genes_url)
    panel_response.raise_for_status()
    panel_json = panel_response.json()

    for gene in panel_json["genes"]:

        # only retain green genes
        if gene["confidence_level"] != "3" or gene["entity_type"] != "gene":
            continue

        ensg = None
        symbol = gene.get("entity_name")

        # take the PanelApp MOI and simplify
        moi = get_simple_moi(gene.get("mode_of_inheritance", None))

        # for some reason the build is capitalised oddly in panelapp
        # at least one entry doesn't have an ENSG annotation
        for build, content in gene["gene_data"]["ensembl_genes"].items():
            if build.lower() == reference_genome.lower():
                # the ensembl version may alter over time, but will be singular
                ensg = content[list(content.keys())[0]]["ensembl_id"]

        # this appears to be missing in latest panel version
        if symbol == 'RNU12' and ensg is None:
            ensg = 'ENSG00000276027'

        # save the entity into the final dictionary
        gene_dict[ensg] = {
            'symbol': symbol,
            'moi': moi,
        }

    return gene_dict


def get_panel_changes(
    previous_version: str, panel_id: str, latest_content: Dict[str, Dict[str, str]]
) -> Dict[str, Dict[str, str]]:
    """
    take the latest panel content, and compare with a previous version
    https://panelapp.agha.umccr.org/api/v1/panels/137/?version=0.10952

    What do I want to get back?
    - entirely new entries
    - changed MOI

    :param previous_version:
    :param panel_id:
    :param latest_content:
    :return:
    """

    updated_content = {'new': [], 'changed': {}}

    # get the full content for the specified panel version
    previous_content = get_panel_green(panel_id=panel_id, version=previous_version)

    # iterate over the latest content
    for key, value in latest_content.items():

        # if the gene wasn't present before, take it in full
        if key not in previous_content:
            updated_content['new'].append(key)

        # otherwise check if the MOI has changed
        else:
            prev_moi = previous_content.get(key).get('moi')
            latest_moi = value.get('moi')

            # if so, store the old and new MOI
            if prev_moi != latest_moi:
                updated_content['changed'][key] = (prev_moi, latest_moi)
    return updated_content


@click.command()
@click.option('--id', 'panel_id', default='137', help='ID to use in panelapp')
@click.option('--out', 'out_path', help='path to write resulting JSON to')
@click.option(
    '--date',
    'since',
    help='identify panel differences between this date and now (YYYY-MM-DD)',
)
def main(panel_id: str, out_path: str, since: Optional[str] = None):
    """
    takes a panel ID and a date
    finds all current panel data from the API
    uses activities endpoint to get highest panel version up to the given date
    retrieves panel data at that point
    compares, and records all panel differences
        - new genes
        - altered MOI
    :param panel_id:
    :param out_path: path to write a JSON object out to
    :param since: string to parse as a Datetime
    :return:
    """

    latest_green = get_panel_green(panel_id=panel_id)
    panel_dict = {'panel_data': latest_green}

    if since is not None:
        since = datetime.strptime(since, "%Y-%m-%d")
        if since > datetime.today():
            raise ValueError(f'The specified date {since} cannot be in the future')
        early_version = get_previous_version(panel_id=panel_id, since=since)
        logging.info('Previous panel version: %s', early_version)
        panel_dict['changes'] = get_panel_changes(
            previous_version=early_version,
            panel_id=panel_id,
            latest_content=latest_green,
        )

    else:
        panel_dict['changes'] = {}

    with open(out_path, 'w', encoding='utf-8') as handle:
        json.dump(panel_dict, handle, indent=True, default=str)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=E1120
