#!/usr/bin/env python3  # noqa: EXE001

"""
This scripts extracts the HPO terms for all individuals in a dataset, as
well as the pedigree and internal SG : external sample ID mapping. These are
saved to files on disk before being zipped and uploaded into the release
bucket of the dataset, in a directory containing today's date.
The latest dataset-vcf type analyses for the dataset are also uploaded.
"""

import csv
import json
import logging
import os
import subprocess
import sys
from datetime import datetime, timezone
from zipfile import ZipFile

import click
from google.cloud import storage
from metamist.apis import ParticipantApi
from metamist.graphql import query, gql

PAPI = ParticipantApi()

TODAY = datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')


def get_individual_metadata(dataset: str):
    """Returns all rows of seqr individual metadata for a dataset"""
    # TODO: replace this with GraphQL once seqr metadata fields are exposed
    return PAPI.get_individual_metadata_for_seqr(dataset).get('rows')


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


def get_pedigrees(dataset: str):
    """Get all pedigree rows for a dataset"""
    _query = gql(
        """
        query DatasetPedigree($datasetName: String!) {
            project(name: $datasetName) {
                pedigree
            }
        }
        """
    )

    return query(_query, {'datasetName': dataset})['project']['pedigree']


def get_sg_id_to_family_guid_map(dataset: str):
    """Reads the json files in the bucket and returns a mapping of SG ID to family GUID"""
    exome_family_guid_map = {}
    with open(
        f'gs://cpg-{dataset}-main-upload/seqr_metadata/udn-aus_exome_seqr_processed.json',
    ) as f:
        exome_family_guid_map = json.load(f)
    
    genome_family_guid_map = {}
    with open(
        f'gs://cpg-{dataset}-main-upload/seqr_metadata/udn-aus_genome_seqr_processed.json',
    ) as f:
        genome_family_guid_map = json.load(f)
        
    return {**exome_family_guid_map, **genome_family_guid_map}


def get_family_guid_map(pedigrees: str, sg_participant_map: dict[str, str], sg_id_family_guid_map: dict[str, str]):
    """Returns a mapping of family ID to GUID"""
    family_guid_map = {}
    for row in pedigrees:
        individual_id = row['individual_id']
        sg_id = sg_participant_map.get(individual_id)
        family_guid = sg_id_family_guid_map.get(sg_id)
        family_guid_map[row['family_id']] = family_guid

    return family_guid_map


def get_participant_sg_map(dataset: str):
    """Returns the internal SG id : external participant id mapping for a dataset."""
    _query = gql(
        """
        query DatasetSampleMap($datasetName: String!) {
            project(name: $datasetName) {
                participants {
                    externalId
                    samples {
                        id
                        externalId
                        sequencingGroups {
                            id
                        }
                    }
                }
            }
        }
        """
    )

    participants = query(_query, {'datasetName': dataset})['project']['participants']

    sg_id_participant_ext_id_map = {}
    for participant in participants:
        participant_ext_id = participant['externalId']
        participant_samples = participant['samples']
        for sample in participant_samples:
            if not sample['sequencingGroups']:
                continue
            sg_id = sample['sequencingGroups'][0]['id']
            sg_id_participant_ext_id_map[sg_id] = participant_ext_id

    return sg_id_participant_ext_id_map


def write_outputs(
    dataset: str,
    individual_hpo_terms: dict[str, list],
    pedigrees: list[dict],
    sg_partitipant_map: dict[str, str],
    family_guid_map: dict[str, str],
    output_path: str,
):
    """Writes the two HPO terms files, pedigree file, and sample map file to a zip."""

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
        json.dump(sg_partitipant_map, f, indent=4, sort_keys=True)
    logging.info(f'Wrote SG ID : Participant external ID map json to {output_path}.')

    # Family GUID map
    with open(f'{output_path}/family_guid_map.json', 'w') as f:
        json.dump(family_guid_map, f, indent=4, sort_keys=True)
    logging.info(f'Wrote family GUID map json to {output_path}.')
    
    # Zip everything
    with ZipFile(f'{dataset}_metadata.zip', 'w') as z:
        z.write(f'{output_path}/hpo_terms.tsv')
        z.write(f'{output_path}/hpo_terms.json')
        z.write(f'{output_path}/pedigree.fam')
        z.write(f'{output_path}/external_translation.json')
        z.write(f'{output_path}/family_guid_map.json')


def upload_metadata_to_release(dataset: str, billing_project: str | None):
    """Uploads the compressed metadata.zip into a directory with today's date in the release bucket"""

    zip_upload_path = os.path.join(TODAY, f'{dataset}_metadata.zip')

    release_bucket = f'cpg-{dataset}-release'

    if billing_project:
        client = storage.Client(project=billing_project)
    else:
        client = storage.Client()

    bucket = client.bucket(release_bucket, user_project=billing_project)

    zip_blob = bucket.blob(zip_upload_path)

    zip_blob.upload_from_filename(f'{dataset}_metadata.zip')

    logging.info(
        f'Uploaded metadata.tar.gz to gs://{os.path.join(release_bucket, zip_upload_path,)}',
    )


def copy_vcf_to_release(dataset: str, billing_project: str | None):
    """Copies the vcf created by the seqr loader to the metadata directory in the release bucket"""
    _query = """
             query AnalysisVCF($datasetName: String!) {
                project(name: $datasetName) {
                    analyses(type: {eq: "CUSTOM"}, status: {eq: COMPLETED}) {
                        id
                        meta
                        output
                        timestampCompleted
                        status
                    }
                }
             }
             """

    analyses = query(_query, {'datasetName': dataset})['project']['analyses']

    # Filter the analyses down to those that have the meta field 'type': 'dataset-vcf'
    vcf_analyses = [
        analysis
        for analysis in analyses
        if ('type', 'dataset-vcf') in analysis['meta'].items()
    ]

    if not vcf_analyses:
        raise RuntimeError(f'{dataset}: No completed dataset-VCF analyses found.')

    vcf_paths = []
    vcf_file_renames = {}

    # Find the latest dataset-vcf analysis based on the timestamp - for both exome and genome
    exome_vcf_analyses = [
        analysis
        for analysis in vcf_analyses
        if analysis['meta'].get('sequencing_type') == 'exome'
    ]

    if exome_vcf_analyses:
        exome_vcf_analyses = sorted(
            exome_vcf_analyses,
            key=lambda a: datetime.strptime(
                a['timestampCompleted'],
                '%Y-%m-%dT%H:%M:%S',
            ).astimezone(),
        )
        latest_exome_analysis = exome_vcf_analyses[-1]

        latest_exome_analysis_date = datetime.strftime(
            datetime.strptime(
                latest_exome_analysis['timestampCompleted'],
                '%Y-%m-%dT%H:%M:%S',
            )
            .astimezone()
            .date(),
            '%Y-%m-%d',
        )

        vcf_file_renames[
            latest_exome_analysis['output']
        ] = f'{latest_exome_analysis_date}_{dataset}_exomes.vcf.bgz'
        vcf_file_renames[
            latest_exome_analysis['output'] + '.tbi'
        ] = f'{latest_exome_analysis_date}_{dataset}_exomes.vcf.bgz.tbi'

        vcf_paths.extend(
            [latest_exome_analysis['output'], latest_exome_analysis['output'] + '.tbi'],
        )
    else:
        logging.info(f'{dataset}: No completed exome VCF analyses found.')

    genome_vcf_analyses = [
        analysis
        for analysis in vcf_analyses
        if analysis['meta'].get('sequencing_type') == 'genome'
    ]
    if genome_vcf_analyses:
        genome_vcf_analyses = sorted(
            genome_vcf_analyses,
            key=lambda a: datetime.strptime(
                a['timestampCompleted'],
                '%Y-%m-%dT%H:%M:%S',
            ).astimezone(),
        )
        latest_genome_analysis = genome_vcf_analyses[-1]

        latest_genome_analysis_date = datetime.strftime(
            datetime.strptime(
                latest_genome_analysis['timestampCompleted'],
                '%Y-%m-%dT%H:%M:%S',
            )
            .astimezone()
            .date(),
            '%Y-%m-%d',
        )

        vcf_file_renames[
            latest_genome_analysis['output']
        ] = f'{latest_genome_analysis_date}_{dataset}_genomes.vcf.bgz'
        vcf_file_renames[
            latest_genome_analysis['output'] + '.tbi'
        ] = f'{latest_genome_analysis_date}_{dataset}_genomes.vcf.bgz.tbi'

        vcf_paths.extend(
            [
                latest_genome_analysis['output'],
                latest_genome_analysis['output'] + '.tbi',
            ],
        )
    else:
        logging.info(f'{dataset}: No completed genome VCF analyses found.')

    # Save the paths to the .vcf.bgz and .vcf.bgz.tbi files and upload them to the release bucket
    if not billing_project:
        billing_project = dataset
    release_path = f'gs://cpg-{dataset}-release/{TODAY}/'
    for vcf_file_path in vcf_paths:
        release_file_path = os.path.join(release_path, vcf_file_renames[vcf_file_path])
        subprocess.run(
            [  # noqa: S603, S607
                'gcloud',
                'storage',
                '--billing-project',
                billing_project,
                'cp',
                vcf_file_path,
                release_file_path,
            ],
            check=True,
        )
    logging.info(f'Copied {vcf_paths} into {release_path}')


@click.command()
@click.option('--dataset')
@click.option('--billing-project', default=None)
@click.option('--metadata-only', is_flag=True)
@click.option('--dry-run', is_flag=True)
def main(dataset: str, billing_project: str | None, metadata_only: bool, dry_run: bool):
    """Creates the metadata files and saves them to the output path"""

    output_path = f'{dataset}_metadata'
    if not os.path.exists(f'./{output_path}'):
        os.makedirs(f'./{output_path}')

    individual_metadata = get_individual_metadata(dataset)
    individual_hpo_terms = get_hpo_terms(individual_metadata)

    pedigrees = get_pedigrees(dataset)

    sg_participant_map = get_participant_sg_map(dataset)
    
    family_guid_map = get_family_guid_map(pedigrees, sg_participant_map)

    write_outputs(
        dataset,
        individual_hpo_terms,
        pedigrees,
        sg_participant_map,
        output_path,
    )

    if not dry_run:
        upload_metadata_to_release(dataset, billing_project)
        if not metadata_only:
            copy_vcf_to_release(dataset, billing_project)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
