#!/usr/bin/env python3  # noqa: EXE001

import logging
import os
import subprocess
import sys
from datetime import datetime

import click
from cloudpathlib.anypath import to_anypath
from cpg_utils.config import get_config
from google.cloud import storage

GENOME_PREFIX = 'mt/'
EXOME_PREFIX = 'exome/mt/'
MT_SUFFIXES = [GENOME_PREFIX, EXOME_PREFIX]
SEQ_TYPE_MAP = {GENOME_PREFIX: 'genome', EXOME_PREFIX: 'exome'}


def get_dataset_mt_timestamps(dataset: str):
    """
    Finds the dataset .mt folders and their creation timestamps by
    searching the /mt/ and /exome/mt/ paths for any .mt folders.
    """
    mt_folders_by_seq_type: dict[str, dict[str, datetime]] = {}

    bucket_path = get_config()['storage'][dataset]['default']
    bucket = storage.Client().get_bucket(bucket_path.removeprefix('gs://'))

    for suffix in MT_SUFFIXES:
        mt_timestamps: dict[str, datetime] | None = {}

        mt_root_folder = to_anypath(os.path.join(bucket_path, suffix))
        if not mt_root_folder.exists():
            continue

        # Look inside the gs://cpg-dataset-main/mt/ and /exome/mt/ paths for .mt folders
        mt_folders = list(mt_root_folder.iterdir())
        mt_folders = [
            str(mt_folder).removeprefix(f'{bucket_path}/') for mt_folder in mt_folders
        ]

        if any(not mt_folder.endswith('.mt/') for mt_folder in mt_folders):
            raise ValueError(f'Invalid .mt folder: {mt_folders}')

        for mt_folder in mt_folders:
            blob = bucket.get_blob(mt_folder)
            if not blob:  # If the .mt folder is not a blob, use the _SUCCESS file
                blob = bucket.get_blob(f'{mt_folder}_SUCCESS')
            try:
                mt_timestamps[mt_folder] = blob.time_created
            except AttributeError:
                logging.warning(f'{dataset}:: Could not get timestamp for: {mt_folder}')
                continue

        if mt_timestamps:
            mt_folders_by_seq_type[SEQ_TYPE_MAP[suffix]] = mt_timestamps

    return mt_folders_by_seq_type


def get_mt_folders_to_delete(
    dataset: str,
    mt_timestamps_dict: dict[str, datetime],
):
    """Takes the list of .mt folders and keeps only the most recently created one"""
    dataset_mts = sorted(mt_timestamps_dict.items(), key=lambda x: x[1])
    dataset_mts_to_delete = [mt[0] for mt in dataset_mts[:-1]]

    logging.info(
        f'{dataset}:: Keeping latest: {dataset_mts[-1][0]}, created: {dataset_mts[-1][1].strftime("%Y-%m-%dT%H:%M:%S")}',
    )

    logging.info(
        f'{dataset}:: Found {len(dataset_mts_to_delete)} .mt entries to delete',
    )
    for mt in dataset_mts_to_delete:
        logging.info(mt)
    return dataset_mts_to_delete


def delete_dataset_mts(dataset: str, dataset_mts_to_delete: list[str]):
    """Takes a list of .mt folders and deletes them entirely"""
    subprocess.run(
        [  # noqa: S603, S607
            'gcloud',
            'storage',
            'rm',
            '--recursive',
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

    for dataset in datasets:
        matrixtables_by_seqtype = get_dataset_mt_timestamps(
            dataset=dataset,
        )
        for sequencing_type, matrixtables in matrixtables_by_seqtype.items():
            logging.info(f'{dataset}:: {sequencing_type.upper()}')
            loads_to_delete = get_mt_folders_to_delete(
                dataset=dataset,
                mt_timestamps_dict=matrixtables,
            )

            if not dry_run:
                delete_dataset_mts(
                    dataset=dataset,
                    dataset_mts_to_delete=loads_to_delete,
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
