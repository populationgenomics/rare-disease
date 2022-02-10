"""
methods for pulling all genes, ENSGs, MOIs,
and activity dates relevant to a gene panel

for each gene on a given panel, store:
- ensg
- symbol
- moi (simplified)
- date of 'Green' status
"""

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


def get_gene_timestamps(panel: Union[str, int]) -> Dict[str, datetime]:
    """
    query for the activities associated with a panel
    associate each gene symbol with the most recent date
    the event date must be the assignment of 'Green'
    :param panel: panel ID to use
    :return:
    """

    # this is a list of dictionaries
    default_date = datetime(1970, 1, 1)
    activity_list = get_json_response(url=ACTIVITIES.format(panel_id=panel))
    date_dict = {}
    for entry in activity_list:

        # only interested in Green gene assignment
        # not worth sophisticated text parsing, just use lowercase
        # "THUMPD1 was set to GREEN", "Rated green"
        if (
            entry.get('entity_type') != 'gene'
            or 'green' not in entry.get('text').lower()
        ):
            continue

        # ensg values not present from this endpoint
        symbol = entry.get('entity_name')

        # uses the django datestamp format
        entry_date = datetime.strptime(entry.get('created'), "%Y-%m-%dT%H:%M:%S.%fz")

        # if the entry date of this activity is more recent, replace/update
        if date_dict.get(symbol, default_date) < entry_date:
            date_dict[symbol] = entry_date

    return date_dict


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
    date_dict: Dict[str, datetime],
    reference_genome: Optional[str] = "GRch38",
    panel_number: str = '137',
) -> Dict[str, Dict[str, Union[datetime, str]]]:
    """
    Takes a panel number, and pulls all details from PanelApp
    :param date_dict: a lookup of gene symbol to entry date onto panel
    :param reference_genome: GRch37 or GRch38
    :param panel_number: defaults to the PanelAppAU Mendeliome
    :return:
    """

    gene_dict = {}

    panel_app_genes_url = PANEL_CONTENT.format(panel_id=panel_number)

    panel_response = requests.get(panel_app_genes_url)
    panel_response.raise_for_status()
    panel_json = panel_response.json()

    # if we don't have a date from activities endpoint, use today
    default_date = datetime.today()

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

        # datetime isn't easily serializable, so save a String repr
        gene_dict[ensg] = {
            'symbol': symbol,
            'ensg': ensg,
            'moi': moi,
            'since': date_dict.get(symbol, default_date).strftime('%Y-%m-%d'),
        }

    return gene_dict


@click.command()
@click.option('--id', 'panel_id', default='137', help='ID to use in panelapp')
@click.option('--out', 'out_path', help='path to write resulting JSON to')
def main(panel_id: str, out_path: str):
    """
    takes a panel ID, and yields a collection of panel data objects
    :param panel_id:
    :param out_path: path to write a JSON object out to
    :return:
    """
    gene_dates = get_gene_timestamps(panel=panel_id)
    panel_genes = get_panel_green(date_dict=gene_dates, panel_number=panel_id)
    with open(out_path, 'w', encoding='utf-8') as handle:
        json.dump(panel_genes, handle, indent=True, default=str)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
