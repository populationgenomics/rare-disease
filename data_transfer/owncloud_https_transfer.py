#!/usr/bin/env python3

"""
Transfer data files from MCRI owncloud cURLs to a dataset's GCP main-upload bucket.
Also works for any other downloads that can be done via "Copy as CURL" command from a browser.

HOWTO:
1. Go to the owncloud page (or other source, e.g. BluePrintGenetics) and right-click > inspect.
2. Open the network tab, then click the download links for the files you want to download.
3. Cancel the download that has started in your broswer.
4. Right-click on the download request in the network tab, and select "Copy as cURL".
5. Paste the cURL command into a file, one per line.
6. Remove all linebreaks and '\ ' characters from the cURL command so each line is a singluar
   command, containing all the necessary information, e.g. cookies, referer, etc.
7. Save the file and run this script with the file path as the argument.
"""

import os

import click
import hailtop.batch as hb
from cloudpathlib import AnyPath
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    dataset_path,
    remote_tmpdir,
)


@click.command('Transfer data from owncloud cURLs')
@click.option('--owncloud-curl-file-path')
def main(owncloud_curl_file_path: str):
    """
    Given a list of cURL commands, download the files and upload them to GCS.
    GCP suffix in target GCP bucket is defined using analysis-runner's --output
    """

    env_config = get_config()
    cpg_driver_image = env_config['workflow']['driver_image']
    billing_project = env_config['hail']['billing_project']
    dataset = env_config['workflow']['dataset']
    output_prefix = env_config['workflow']['output_prefix']
    assert all({billing_project, cpg_driver_image, dataset, output_prefix})

    with AnyPath(owncloud_curl_file_path).open() as file:
        owncloud_curls = [line.strip() for line in file.readlines() if line.strip()]

    incorrect_urls = [url for url in owncloud_curls if not url.startswith('\'https://')]
    if incorrect_urls:
        raise Exception(f'Incorrect URLs: {incorrect_urls}')

    sb = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(f'transfer {dataset}', backend=sb, default_image=cpg_driver_image)
    output_path = dataset_path(output_prefix, 'upload')

    # may as well batch them to reduce the number of VMs
    for idx, curl in enumerate(owncloud_curls):
        if curl.startswith('curl '):
            curl.removeprefix('curl ')
        url = curl.split(' ')[0]
        try:
            filename = os.path.basename(url).split('&files=')[1].removesuffix("'")
        except IndexError:
            filename = url.split('?')[0].split('/')[-1]
        if '&downloadStartSecret' in filename:
            filename = filename.split('&downloadStartSecret')[0]
        if not filename:
            filename = f'file_{idx}.tar'
        j = batch.new_job(f'URL {idx} ({filename})')
        authenticate_cloud_credentials_in_job(job=j)
        # catch errors during the cURL
        j.command('set -euxo pipefail')
        j.command(f'curl -L {curl} | gsutil cp - {os.path.join(output_path, filename)}')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
    # flake8: noqa
