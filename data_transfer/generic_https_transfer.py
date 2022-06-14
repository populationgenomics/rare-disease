#!/usr/bin/env python3
"""
Transfer datasets from presigned URLs to a dataset's GCP main-upload bucket.
"""
import os
from shlex import quote

import click
import hailtop.batch as hb
from cloudpathlib import AnyPath
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    dataset_path,
    remote_tmpdir,
)


@click.command("Transfer_datasets from signed URLs")
@click.option("--presigned-url-file-path")
def main(presigned_url_file_path: str):
    """
    Given a list of presigned URLs, download the files and upload them to GCS.
    GCP suffix in target GCP bucket is defined using analysis-runner's --output
    """

    env_config = get_config()
    cpg_driver_image = env_config["workflow"]["driver_image"]
    billing_project = env_config["hail"]["billing_project"]
    dataset = env_config["workflow"]["dataset"]
    output_prefix = env_config["workflow"]["output_prefix"]
    assert all({billing_project, cpg_driver_image, dataset, output_prefix})

    with open(AnyPath(presigned_url_file_path)) as file:
        presigned_urls = [l.strip() for l in file.readlines() if l.strip()]

    incorrect_urls = [url for url in presigned_urls if not url.startswith("https://")]
    if incorrect_urls:
        raise Exception(f"Incorrect URLs: {incorrect_urls}")

    sb = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(f"transfer {dataset}", backend=sb, default_image=cpg_driver_image)

    output_path = dataset_path(output_prefix, "upload")

    # may as well batch them to reduce the number of VMs
    for idx, url in enumerate(presigned_urls):

        filename = os.path.basename(url).split("?")[0]
        j = batch.new_job(f"URL {idx} ({filename})")
        quoted_url = quote(url)
        authenticate_cloud_credentials_in_job(job=j)
        # catch errors during the cURL
        j.command("set -euxo pipefail")
        j.command(
            f"curl -L {quoted_url} | gsutil cp - {os.path.join(output_path, filename)}"
        )

    batch.run(wait=False)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
