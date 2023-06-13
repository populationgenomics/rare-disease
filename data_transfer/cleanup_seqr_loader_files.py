#!/usr/bin/env python3  # noqa: EXE001

import logging
import subprocess
import sys
from datetime import datetime, timedelta, timezone

import click
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
        if (seqr_dir := blob.name.removeprefix(prefix).split('/')[0]) not in seqr_load_directories:
            seqr_load_directories.add(seqr_dir)

        if blob.name.endswith('full.vcf.gz'):
            load_id = blob.name.removesuffix(VCF_SUFFIX)
            seqr_loads[load_id] = blob.time_created

    logging.info(f'Found {len(seqr_load_directories)} seqr load directories')
    logging.info(f'Found {len(seqr_loads)} seqr loads with full VCF')

    return seqr_loads, seqr_load_directories


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

    seqr_loads_to_delete = [
        f'gs://cpg-seqr-main/{load}' for load in seqr_loads_to_delete
    ]
    subprocess.run(
        [  # noqa: S603, S607
            'gsutil',
            'rm',
            '-f',
            '-r',
            *seqr_loads_to_delete,
        ],
        check=True,
    )
    logging.info(f'Deleted: {seqr_loads_to_delete}')


@click.command()
@click.option('--dry-run', is_flag=True)
def main(dry_run: bool):
    """Run the seqr_load analyser and deleter for genomes and exomes"""

    logging.info('GENOME')
    genome_loads, extra_genome_folders = get_seqr_loads(
        bucket=BUCKET,
        prefix=GENOME_PREFIX,
    )

    for folder in extra_genome_folders:
        logging.info(f'Genome folder {folder} has no VCF')

    genome_loads_to_delete = get_seqr_loads_to_delete(seqr_loads_dict=genome_loads)

    if not dry_run:
        delete_seqr_load_files(genome_loads_to_delete)

    logging.info('EXOME')
    exome_loads, extra_exome_folders = get_seqr_loads(
        bucket=BUCKET,
        prefix=EXOME_PREFIX,
    )

    for folder in extra_exome_folders:
        logging.info(f'Exome folder {folder} has no VCF')

    exome_loads_to_delete = get_seqr_loads_to_delete(seqr_loads_dict=exome_loads)

    if not dry_run:
        delete_seqr_load_files(exome_loads_to_delete)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    # pylint: disable=no-value-for-parameter
    main()
