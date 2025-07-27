#!/usr/bin/env python3  # noqa: EXE001
"""
  !!    THIS SCRIPT CAN DELETE PRODUCTION DATA. USE WITH CAUTION.     !!

  - Seqr bucket matrixtable cleanup script -
This script is used to delete old matrixtables from the seqr-main bucket.
These tables are created with the AnnotateCohort step of the seqr_loader pipeline, for both exomes and genomes.
Any matrixtables that are older than 30 days will be deleted, unless they are the most recent one.

Once AnnotateCohort has its own Metamist analysis type, this script can be replaced with one that
relies on Metamist analysis completion timestamps instead of the matrixtable creation timestamps.
"""
import logging
import os
import sys
from datetime import datetime, timedelta, timezone

import click
from cloudpathlib import CloudPath
from cpg_utils.config import get_config
from google.cloud import storage

TODAY = datetime.now(tz=timezone.utc)
CLIENT = storage.Client()
BUCKET = CLIENT.get_bucket('cpg-seqr-main')

GENOME_PREFIX = 'seqr_loader/'
EXOME_PREFIX = 'exome/seqr_loader/'
VCF_SUFFIX = '/JointGenotyping/full.vcf.gz'


def get_seqr_loads(bucket, prefix: str):  # noqa: ANN001
    """Finds the seqr_loader directories and VCF create times in the seqr-main bucket"""

    seqr_load_directories = set()
    seqr_loads = {}
    for blob in bucket.list_blobs(prefix=prefix):
        if '.mt' in blob.name or '.ht' in blob.name:
            continue

        # Some directories have no VCF, but we still want them
        if (
            seqr_dir := blob.name.removeprefix(prefix).split('/')[0]
        ) not in seqr_load_directories:
            seqr_load_directories.add(seqr_dir)

        if blob.name.endswith('full.vcf.gz'):
            load_id = blob.name.removesuffix(VCF_SUFFIX)
            seqr_loads[load_id] = blob.time_created

    logging.info(f'Found {len(seqr_load_directories)} seqr load directories')
    logging.info(f'Found {len(seqr_loads)} seqr loads with full VCF')

    load_folders = {path.removeprefix(prefix) for path in seqr_loads}
    extra_folders = seqr_load_directories.difference(load_folders)

    return seqr_loads, extra_folders


def get_seqr_loads_to_delete(seqr_loads_dict: dict):
    """Takes the list of VCFs from the seqr loads folder, keeping the latest VCF and any younger than 30 days"""

    # Keep the latest seqr load regardless of age
    seqr_loads = sorted(seqr_loads_dict.items(), key=lambda x: x[1])
    latest_seqr_load = seqr_loads[-1]
    seqr_loads.remove(latest_seqr_load)

    logging.info(
        f'Keeping latest seqr load: {latest_seqr_load[0]}, created: {latest_seqr_load[1].strftime("%Y-%m-%d")}',
    )

    seqr_loads_to_delete = []
    for seqr_load, time_created in seqr_loads:
        if time_created < TODAY - timedelta(days=30):
            seqr_loads_to_delete.append(seqr_load)

    logging.info(f'Found {len(seqr_loads_to_delete)} seqr loads to delete')
    return seqr_loads_to_delete


def delete_seqr_load_files(seqr_loads_to_delete: list[str]):
    """Takes a list of seqr_loader folders and deletes them entirely"""
    config = get_config()
    bucket_path = config['storage']['seqr']['default']
    seqr_loads_to_delete = [
        os.path.join(bucket_path, load) for load in seqr_loads_to_delete
    ]
    for seqr_load in seqr_loads_to_delete:
        seqr_load_path = CloudPath(seqr_load)
        seqr_load_path.rmtree()
        logging.info(f'Deleted: {seqr_load}')


@click.command()
@click.option('--dry-run', is_flag=True)
def main(dry_run: bool):
    """Run the seqr_load analyser and deleter for genomes and exomes"""
    matrixtables_to_delete = []

    for prefix in [GENOME_PREFIX, EXOME_PREFIX]:
        logging.info(f'Prefix: {prefix}')
        seqr_loads, extra_folders = get_seqr_loads(
            bucket=BUCKET,
            prefix=prefix,
        )

        for folder in extra_folders:
            logging.info(f'Folder {folder} has no VCF')

        matrixtables_to_delete.extend(
            get_seqr_loads_to_delete(seqr_loads_dict=seqr_loads),
        )

    if dry_run:
        logging.info('<< Dry run, not deleting anything >>')
        logging.info(
            f'Would have deleted {len(matrixtables_to_delete)} matrixtables',
        )
    else:
        delete_seqr_load_files(matrixtables_to_delete)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
