#!/usr/bin/env python3
"""
Transfer datasets from presigned URLs to a dataset's GCP main-upload bucket.
"""
from typing import Optional
import os
from shlex import quote

import click
import hailtop.batch as hb
from cloudpathlib import AnyPath
from cpg_utils.hail import remote_tmpdir
from analysis_runner.constants import GCLOUD_ACTIVATE_AUTH


DRIVER_IMAGE = os.getenv("CPG_DRIVER_IMAGE")
DATASET = os.getenv("CPG_DATASET")

assert DRIVER_IMAGE and DATASET


@click.command("Transfer_datasets from signed URLs")
@click.option("--presigned-url-file-path")
@click.option("--subfolder", type=str)
def main(
    presigned_url_file_path: str,
    subfolder: Optional[str] = None,
):
    """
    Given a list of presigned URLs, download the files and upload them to GCS.
    """

    with open(AnyPath(presigned_url_file_path)) as file:
        presigned_urls = [l.strip() for l in file.readlines() if l.strip()]

    incorrect_urls = [url for url in presigned_urls if not url.startswith("https://")]
    if incorrect_urls:
        raise Exception(f"Incorrect URLs: {incorrect_urls}")

    sb = hb.ServiceBackend(
        billing_project=os.getenv("HAIL_BILLING_PROJECT"),
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(f"transfer {DATASET}", backend=sb, default_image=DRIVER_IMAGE)

    output_path = f"gs://cpg-{DATASET}-main-upload"
    if subfolder:
        output_path = os.path.join(output_path, subfolder)

    # may as well batch them to reduce the number of VMs
    for idx, url in enumerate(presigned_urls):

        filename = os.path.basename(url).split("?")[0]
        j = batch.new_job(f"URL {idx} ({filename})")
        quoted_url = quote(url)
        j.command(GCLOUD_ACTIVATE_AUTH)
        # catch errors during the cURL
        j.command("set -euxo pipefail")
        j.command(
            f"curl -L {quoted_url} | gsutil cp - {os.path.join(output_path, filename)}"
        )

    batch.run(wait=False)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
