#!/usr/bin/env python3  # noqa: EXE001
"""
  !!    THIS SCRIPT CAN DELETE PRODUCTION DATA. USE WITH CAUTION.     !!

  - Bucket matrixtable cleanup script -
This script should be one time use only, and is not intended to be run regularly.
Subsequent cleanups of the dataset matrix tables should rely on Metamist entries.

Find all matrixtables inside each dataset's main bucket.
Searches the /mt/ and /exome/mt/ paths for any .mt folders.
Keeps only the most recently created .mt folder per sequencing type, deletes the rest.
"""
import logging
import os
import sys
from datetime import datetime, timezone

import click
from cloudpathlib import CloudPath
from cloudpathlib.anypath import to_anypath
from cpg_utils.config import get_config
from google.cloud import storage

TODAY = datetime.now(tz=timezone.utc)

GENOME_PREFIX = 'mt/'
EXOME_PREFIX = 'exome/mt/'
MT_PREFIXES = [GENOME_PREFIX, EXOME_PREFIX]
SEQ_TYPE_MAP = {GENOME_PREFIX: 'genome', EXOME_PREFIX: 'exome'}


def get_dataset_mt_timestamps(dataset: str):
    """
    Finds the dataset .mt folders and their creation timestamps by
    searching the /mt/ and /exome/mt/ paths for any .mt folders.
    """
    mt_folders_by_seq_type: dict[str, dict[str, datetime]] = {}

    bucket_path = get_config()['storage'][dataset]['default']
    bucket = storage.Client().get_bucket(bucket_path.removeprefix('gs://'))

    for prefix in MT_PREFIXES:
        mt_timestamps: dict[str, datetime] | None = {}

        mt_root_folder = to_anypath(os.path.join(bucket_path, prefix))
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
                mt_timestamps[
                    os.path.join(bucket_path, prefix, mt_folder)
                ] = blob.time_created
            except AttributeError:
                logging.warning(f'{dataset}:: Could not get timestamp for: {mt_folder}')
                continue

        if mt_timestamps:
            mt_folders_by_seq_type[SEQ_TYPE_MAP[prefix]] = mt_timestamps

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


def write_dataset_mts_to_delete_file(dataset_mts_to_delete: list[str]):
    """Takes the list of matrixtables to delete and writes them to a file in the seqr bucket"""
    bucket_path = get_config()['storage']['seqr']['default']
    bucket = storage.Client().get_bucket(bucket_path.removeprefix('gs://'))

    blob = bucket.blob(f'mt/matrixtables_to_delete_{TODAY}.txt')
    blob.upload_from_string('\n'.join(dataset_mts_to_delete))

    logging.info(
        f'Wrote {len(dataset_mts_to_delete)} matrixtables to delete to: {blob.path}',
    )


@click.command()
@click.option('--dry-run', is_flag=True)
@click.option(
    '--datasets',
    '-d',
    multiple=True,
)
def main(dry_run: bool, datasets: list[str]):
    """Find .mt folders in dataset buckets for genomes and exomes"""
    matrixtables_to_delete = []
    for dataset in datasets:
        matrixtables_by_seqtype = get_dataset_mt_timestamps(
            dataset=dataset,
        )
        for sequencing_type, matrixtables in matrixtables_by_seqtype.items():
            logging.info(f'{dataset}:: {sequencing_type.upper()}')
            matrixtables_to_delete.extend(
                get_mt_folders_to_delete(
                    dataset=dataset,
                    mt_timestamps_dict=matrixtables,
                ),
            )

    write_dataset_mts_to_delete_file(matrixtables_to_delete)
    if dry_run:
        logging.info('<< Dry run, not deleting anything >>')
        logging.info(
            f'Would have deleted {len(matrixtables_to_delete)} .mt analysis outputs from {len(datasets)} datasets',
        )
    else:
        for matrixtable in matrixtables_to_delete:
            logging.info(f'Deleting {matrixtable}')
            mt = CloudPath(matrixtable)
            mt.rmtree()

        logging.info(
            f'Deleted {len(matrixtables_to_delete)} .mt analysis outputs from {len(datasets)} datasets',
        )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
