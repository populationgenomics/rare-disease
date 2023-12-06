"""
script to crawl around in panelapp and extract all the pmids
"""

import json
from time import sleep

import requests

PANELAPP_URL = 'https://panelapp.agha.umccr.org/api/v1/panels'


def get_json_response(
    url: str,
    max_retries: int = 4,
    base_delay: int = 1,
    max_delay: int = 32,
):
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    Args:
        url (str): URL to retrieve JSON format data from
        max_retries (int): maximum number of retries
        base_delay (int): initial delay between retries
        max_delay (int): maximum delay between retries

    Returns:
        the JSON response from the endpoint
    """

    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(
                url,
                headers={'Accept': 'application/json'},
                timeout=60,
            )
            response.raise_for_status()  # Raise an exception for bad responses (4xx and 5xx)
            return response.json()
        except (requests.RequestException, TimeoutError):
            retries += 1
            if retries < max_retries:
                delay = min(base_delay * 2**retries, max_delay)
                sleep(delay)

    raise TimeoutError('Max retries reached. Request failed.')


def find_all_panels() -> dict[int, dict[str, str | int | float | dict]]:
    """
    iterate over all panels, get IDs, names, descriptions, versions
    """
    panel_dict = {}
    panel_url = PANELAPP_URL
    while panel_url:
        panel_data = get_json_response(panel_url)
        for panel in panel_data['results']:
            panel_dict[panel['id']] = {
                'name': panel['name'],
                'description': panel['disease_group'],
                'version': panel['version'],
                'disorders': panel['relevant_disorders'],
                'genes': {},
            }
        panel_url = panel_data['next']
    return panel_dict


def get_ensembl_gene_id(gene_data: dict) -> str | None:
    """
    ronseal

    "ensembl_genes": {
          "GRch37": ...,
          "GRch38": {
            "90": {
              "location": "1:161766294-161964070",
              "ensembl_id": "ENSG00000118217"
            }
          }
        },

    Args:
        gene_data (dict) ensembl gene data
    """
    if (
        (ensembl := gene_data['ensembl_genes'])
        and 'GRch38' in ensembl
        and '90' in ensembl['GRch38']
    ):
        return gene_data['ensembl_genes']['GRch38']['90']['ensembl_id']
    return None


def get_panel_gene_pmids(panel_id: int) -> dict[str, dict]:
    """
    return a dictionary of genes and their pmids, evidence level
    Args:
        panel_id (int): the panel ID to hit
    """
    gene_dict = {}
    panel_data = get_json_response(f'{PANELAPP_URL}/{panel_id}')
    for gene in panel_data['genes']:
        gene_data = gene['gene_data']

        # shove this in at the top level
        gene_data['ensembl_id'] = get_ensembl_gene_id(gene_data)

        gene_dict[gene_data['gene_symbol']] = {
            key: gene_data[key]
            for key in [
                'alias',
                'biotype',
                'gene_symbol',
                'ensembl_id',
                'hgnc_id',
                'omim_gene',
            ]
        } | {
            key: gene[key]
            for key in [
                'confidence_level',
                'entity_name',
                'entity_type',
                'mode_of_inheritance',
                'phenotypes',
                'publications',
            ]
        }

    return gene_dict


if __name__ == '__main__':
    # find all the panels!
    panel_dict = find_all_panels()

    # iterate over all the panels, get the genes
    for panel, content in panel_dict.items():
        content['genes'] = get_panel_gene_pmids(panel)
        sleep(2)

    out_path = 'this_local_filename.json'
    with open(out_path, 'w', encoding='utf-8') as handle:
        json.dump(panel_dict, handle, indent=4)
