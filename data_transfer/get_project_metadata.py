import csv
import json
import logging
import os
import subprocess
import sys
from datetime import datetime, timezone

import click
from sample_metadata.apis import FamilyApi, ParticipantApi, SampleApi
from sample_metadata.graphql import query

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
    _query = """
             query ProjectPedigree($projectName: String!) {
                 project(name: $projectName) {
                     pedigree
                 }
             }
             """

    return query(_query, {'projectName': project}).get('project').get('pedigree')


def get_sample_map(project: str):
    """Returns the internal id : external id sample map for a project"""
    _query = """
            query ProjectSampleMap($projectName: String!) {
                project(name: $projectName) {
                    participants {
                        samples {
                            id
                            externalId
                        }
                    }
                }
            }
            """

    samples_list = (
        query(_query, {'projectName': project}).get('project').get('participants')
    )

    samples = [
        sample
        for samples_dict in samples_list
        for sample in samples_dict.get('samples')
        if samples_dict.get('samples')
    ]

    return {sample.get('id'): sample.get('externalId') for sample in samples if sample}


def write_outputs(
    individual_hpo_terms: dict[str, list],
    pedigrees: list[dict],
    sample_map: dict[str, str],
    output_path: str,
):
    """Writes the HPO terms file, pedigree file, and sample map file to the output path"""

    # HPO terms tsv
    with open(f'{output_path}/hpo_terms.tsv', 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        for individual, hpo_terms in individual_hpo_terms.items():
            row = [individual, *list(hpo_terms)]
            if row:
                writer.writerow(row)
    logging.info(f'Wrote HPO terms tsv to {output_path}.')

    # HPO terms json
    with open(f'{output_path}/hpo_terms.json', 'w') as f:
        json.dump(individual_hpo_terms, f, indent=4, sort_keys=True)
    logging.info(f'Wrote HPO terms json to {output_path}.')

    # Pedigree file
    with open(f'{output_path}/pedigree.fam', 'w', encoding='utf-8') as f:
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
        json.dump(sample_map, f, indent=4, sort_keys=True)
    logging.info(f'Wrote sample ID map json to {output_path}.')

    subprocess.run(
        ['tar', '-zcvf', 'metadata.tar.gz', f'{output_path}'],  # noqa: S603, S607
        check=True,
        stderr=subprocess.DEVNULL,
    )
    logging.info('Compressed metadata into metadata.tar.gz.')


def upload_metadata_to_release(project: str):
    """Uploads the compressed metadata .tar.gz into a directory with today's date in the release bucket"""
    today = datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')
    release_path = f'gs://cpg-{project}-release/{today}/'

    subprocess.run(
        [  # noqa: S603, S607
            'gcloud',
            'storage',
            '--billing-project',
            project,
            'cp',
            './metadata.tar.gz',
            release_path,
        ],
        check=True,
    )
    logging.info(f'Uploaded metadata.tar.gz to {release_path}')


@click.command()
@click.option('--project')
@click.option('--dry-run', is_flag=True)
def main(project: str, dry_run: bool):
    """Creates the metadata files and saves them to the output path"""

    output_path = f'{project}_metadata'
    if not os.path.exists(f'./{output_path}'):
        os.makedirs(f'./{output_path}')

    individual_metadata = get_individual_metadata(project)
    individual_hpo_terms = get_hpo_terms(individual_metadata)

    pedigrees = get_pedigrees(project)

    sample_map = get_sample_map(project)

    write_outputs(individual_hpo_terms, pedigrees, sample_map, output_path)

    if not dry_run:
        upload_metadata_to_release(project)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
