#!/usr/bin/env python3

"""Checks data transfer integrity for GDR cohorts by comparing MD5 checksums."""
import base64
import binascii
import csv
import logging
import os
import sys

import click
from cloudpathlib import AnyPath
from cpg_utils.config import get_config
from google.cloud import storage

try:
    DEFAULT_DATASET = get_config()['workflow']['dataset']
except (KeyError, AssertionError):
    DEFAULT_DATASET = None


def main(
    dataset: str,
    file_prefix: str = '',
    filename_column: str = 'filename',
    checksum_column: str = 'checksum',
    manifest_file_path: str = 'manifest.txt',
    delimiter: str = '\t',
):
    """Main entrypoint."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    upload_bucket_name = f'cpg-{dataset}-main-upload'

    if not manifest_file_path.startswith('gs://'):
        # empty file_prefix gets skipped
        manifest_file_path = os.path.join(
            f'gs://{upload_bucket_name}',
            file_prefix,
            manifest_file_path,
        )

    logging.info(f'Fetching manifest from {manifest_file_path}')
    with AnyPath(manifest_file_path).open() as manifest_file:
        manifest = manifest_file.read()

    storage_client = storage.Client()
    bucket = storage_client.get_bucket(upload_bucket_name)

    any_errors = False
    tsv_reader = csv.DictReader(manifest.splitlines(), delimiter=delimiter)
    matches = 0
    mismatches = 0
    for row in tsv_reader:
        filename = row[filename_column]
        expected_md5 = row[checksum_column]
        # empty file_prefix gets skipped
        check_blob = bucket.get_blob(os.path.join(file_prefix, filename))
        if not check_blob:
            logging.error(f'blob does not exist: {filename}')
            any_errors = True
            continue
        # Read the checksum from the blob. The checksum is base64-encoded.
        actual_md5 = binascii.hexlify(
            base64.urlsafe_b64decode(check_blob.md5_hash),
        ).decode('utf-8')

        if expected_md5 == actual_md5:
            logging.info(f'match: {filename}')
            matches += 1
        else:
            logging.error(f'mismatch: {filename}, {expected_md5=}, {actual_md5=}')
            any_errors = True
            mismatches += 1

    logging.info(f'{matches=}, {mismatches=}')
    if any_errors:
        sys.exit(1)


@click.command()
@click.option(
    '--dataset',
    required=DEFAULT_DATASET is None,
    default=DEFAULT_DATASET,
)
@click.option('--file-prefix', default='')
@click.option(
    '--manifest-file-path',
    required=False,
    default='manifest.txt',
)
@click.option(
    '--filename-column',
    default='filename',
)
@click.option(
    '--checksum-column',
    default='checksum',
)
def from_cli(
    dataset: str,
    file_prefix: str = '',
    filename_column: str = 'filename',
    checksum_column: str = 'checksum',
    manifest_file_path: str = 'manifest.txt',
):
    """Run the script from the command line."""
    return main(
        dataset=dataset,
        file_prefix=file_prefix,
        manifest_file_path=manifest_file_path,
        filename_column=filename_column,
        checksum_column=checksum_column,
    )


if __name__ == '__main__':
    from_cli()  # pylint: disable=no-value-for-parameter
