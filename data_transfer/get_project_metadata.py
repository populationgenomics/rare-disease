#!/usr/bin/env python

import csv
import json
import logging
import sys

import click
from sample_metadata.apis import FamilyApi, ParticipantApi, SampleApi

PAPI = ParticipantApi()
SAPI = SampleApi()
FAPI = FamilyApi()


def get_individual_metadata(project: str):
    """Returns all rows of seqr individual metadata for a project"""
    return PAPI.get_individual_metadata_for_seqr(project).get('rows')


def get_hpo_terms(individual_metadata: list[dict]):
    """Returns a mapping of individual ID to list of hpo terms"""
    indiv_hpo_terms = {}
    for row in individual_metadata:
        indiv_id = row.get('individual_id')
        hpo_terms = ''
        if row.get('hpo_terms_present'):
            hpo_terms = row.get('hpo_terms_present').split(',')

        indiv_hpo_terms[indiv_id] = hpo_terms

    return indiv_hpo_terms


def get_pedigrees(project: str):
    """Get all pedigree rows for a project"""
    return FAPI.get_pedigree(project)


def get_sample_map(project: str):
    """Returns the internal id : external id sample map for a project"""
    return SAPI.get_all_sample_id_map_by_internal(project)


def write_outputs(
    individual_hpo_terms: dict[str, list],
    pedigrees: list[dict],
    sample_map: dict[str, str],
    output_path: str,
):
    """Writes the HPO terms file, pedigree file, and sample map file to the output path"""

    # HPO terms tsv
    with open(f'{output_path}/hpo_terms.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for individual, hpo_terms in individual_hpo_terms.items():
            row = [individual, *list(hpo_terms)]
            if row:
                writer.writerow(row)
    logging.info(f'Wrote HPO terms tsv to {output_path}.')

    # HPO terms json
    with open(f'{output_path}/hpo_terms.json', 'w') as f:
        json.dump(individual_hpo_terms, f, indent=4)
    logging.info(f'Wrote HPO terms json to {output_path}.')

    # Pedigree file
    with open(f'{output_path}/pedigree.fam', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(
            [
                'Family ID',
                'Individual ID',
                'Paternal ID',
                'Maternal ID',
                'Sex',
                'Affected',
            ],
        )
        for ped_row in pedigrees:
            values = list(ped_row.values())
            writer.writerow(values)
    logging.info(f'Wrote pedigree to {output_path}.')

    # Sample map
    with open(f'{output_path}/external_translation.json', 'w') as f:
        json.dump(sample_map, f, indent=4)
    logging.info(f'Wrote sample ID map json to {output_path}.')


@click.command()
@click.option('--project')
@click.option('--output-path')
def main(project: str, output_path: str):
    """Creates the metadata files and saves them to the output path"""
    individual_metadata = get_individual_metadata(project)
    individual_hpo_terms = get_hpo_terms(individual_metadata)

    pedigrees = get_pedigrees(project)
    print(type(pedigrees[0]))

    sample_map = get_sample_map(project)

    write_outputs(individual_hpo_terms, pedigrees, sample_map, output_path)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
