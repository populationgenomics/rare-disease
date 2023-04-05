#!/usr/bin/env python3

"""Checks data transfer integrity for GDR cohorts by comparing MD5 checksums."""

import base64
import binascii
import csv
import logging
import sys
from google.cloud import storage

from cpg_utils.config import get_config

MANIFEST_FILE = 'manifest.txt'
FILENAME_COLUMN = 'filename'
CHECKSUM_COLUMN = 'checksum'


def main():
    """Main entrypoint."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    dataset = get_config()['workflow']['dataset']
    upload_bucket_name = f'cpg-{dataset}-main-upload'

    storage_client = storage.Client()
    bucket = storage_client.get_bucket(upload_bucket_name)

    manifest_blob = bucket.get_blob(MANIFEST_FILE)
    if not manifest_blob:
        logging.error(f'blob does not exist: {MANIFEST_FILE}')
        sys.exit(1)
    manifest = manifest_blob.download_as_text()

    any_errors = False
    tsv_reader = csv.DictReader(manifest.splitlines(), delimiter='\t')
    for row in tsv_reader:
        filename = row[FILENAME_COLUMN]
        expected_md5 = row[CHECKSUM_COLUMN]

        # Read the checksum from the blob. The checksum is base64-encoded.
        check_blob = bucket.get_blob(filename)
        if not check_blob:
            logging.error(f'blob does not exist: {filename}')
            any_errors = True
            continue
        actual_md5 = binascii.hexlify(
            base64.urlsafe_b64decode(check_blob.md5_hash)
        ).decode('utf-8')

        if expected_md5 == actual_md5:
            logging.info(f'match: {filename}')
        else:
            logging.error(f'mismatch: {filename}, {expected_md5=}, {actual_md5=}')
            any_errors = True

    if any_errors:
        sys.exit(1)


if __name__ == '__main__':
    main()
