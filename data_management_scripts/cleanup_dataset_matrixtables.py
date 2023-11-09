#!/usr/bin/env python3  # noqa: EXE001

import logging
import subprocess
import sys
from datetime import datetime

import click
from cpg_utils.config import get_config
from google.cloud import storage

GENOME_PREFIX = 'mt/'
EXOME_PREFIX = 'exome/mt/'
VCF_SUFFIX = '/JointGenotyping/full.vcf.gz'


def get_dataset_mt_timestamps(
    dataset: str,
    bucket: storage.Bucket,
    prefix: str,
):
    """Finds the .mt directories and their create times in the seqr-main bucket"""
    mt_timestamps = {}
    mt_folders = (
        subprocess.run(  # Use subprocess because list_blobs lists everything inside each .mt folder
            args=['gsutil', 'ls', f'gs://{bucket.name}/{prefix}'],
            check=True,
            capture_output=True,
            text=True,
        )
        .stdout.strip()
        .split('\n')
    )

    mt_folders = [
        mt_folder.removeprefix(f'gs://{bucket.name}/') for mt_folder in mt_folders
    ]

    for mt_folder in mt_folders:
        logging.info(f'{dataset}:: Found: {mt_folder}')
        blob = bucket.get_blob(mt_folder)
        if not blob:  # If the .mt folder is not a blob, use the _SUCCESS file
            blob = bucket.get_blob(f'{mt_folder}_SUCCESS')
        mt_timestamps[mt_folder] = blob.time_created

    return mt_timestamps


def get_mt_folders_to_delete(
    dataset: str,
    mt_timestamps_dict: dict[str, datetime],
):
    """Takes the list of .mt directories and keeps only the most recently created one"""
    dataset_mts = sorted(mt_timestamps_dict.items(), key=lambda x: x[1])
    dataset_mts_to_delete = [mt[0] for mt in dataset_mts[:-1]]

    logging.info(
        f'{dataset}:: Keeping latest: "{dataset_mts[-1][0]}", created: {dataset_mts[-1][1].strftime("%Y-%m-%dT%H:%M:%S")}',
    )

    logging.info(
        f'{dataset}:: Found {len(dataset_mts_to_delete)} .mt entries to delete',
    )
    return dataset_mts_to_delete


def delete_dataset_mts(dataset: str, dataset_mts_to_delete: list[str]):
    """Takes a list of .mt folders and deletes them entirely"""
    subprocess.run(
        [  # noqa: S603, S607
            'gsutil',
            '-m',
            'rm',
            '-f',
            '-r',
            *dataset_mts_to_delete,
        ],
        check=True,
    )
    for mt in dataset_mts_to_delete:
        logging.info(f'{dataset}:: Deleted: {mt}')


@click.command()
@click.option('--dry-run', is_flag=True)
@click.option(
    '--datasets',
    '-d',
    multiple=True,
)
def main(dry_run: bool, datasets: list[str]):
    """Run the seqr_load analyser and deleter for genomes and exomes"""
    client = storage.Client()
    config = get_config()

    for dataset in datasets:
        bucket_name = config['storage'][dataset]['default'].removeprefix('gs://')
        bucket = client.get_bucket(bucket_name)

        logging.info(f'{dataset}:: GENOME')
        genome_loads = get_dataset_mt_timestamps(
            dataset=dataset,
            bucket=bucket,
            prefix=GENOME_PREFIX,
        )

        genome_loads_to_delete = get_mt_folders_to_delete(
            dataset=dataset,
            mt_timestamps_dict=genome_loads,
        )

        if not dry_run:
            delete_dataset_mts(
                dataset=dataset,
                dataset_mts_to_delete=genome_loads_to_delete,
            )

        logging.info(f'{dataset}:: EXOME')
        exome_loads = get_dataset_mt_timestamps(
            dataset=dataset,
            bucket=bucket,
            prefix=EXOME_PREFIX,
        )

        exome_loads_to_delete = get_mt_folders_to_delete(
            dataset=dataset,
            mt_timestamps_dict=exome_loads,
        )

        if not dry_run:
            delete_dataset_mts(
                dataset=dataset,
                dataset_mts_to_delete=exome_loads_to_delete,
            )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
