#!/usr/bin/env python3

"""
Use VEP using a dataproc cluster
"""


import os
import hailtop.batch as hb
from analysis_runner import dataproc
import click


@click.command()
@click.option('--script', 'script', help='path to VEP main script')
@click.option('--matrix', required=True, help='Hail matrix table to run VEP on')
@click.option('--output', required=True, help='Hail matrix table to write')
def main(script: str, matrix: str, output: str):
    """
    runs a script inside dataproc to execute VEP
    :param script: str, the path to the VEP main script
    :param matrix: str, the bucket path to a MatrixTable
    :param output: str, the bucket path to a MatrixTable for writing out
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(name='run_vep_in_dataproc_cluster', backend=service_backend)

    job = dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        script=f'{script} --matrix {matrix} --output {output}',
        max_age='8h',
        init=[
            'gs://cpg-reference/hail_dataproc/install_common.sh',
            'gs://cpg-reference/vep/vep-GRCh38-perl5_set.sh',  # install and configure VEP 105
        ],
        job_name='run_vep',
        num_secondary_workers=10,
        num_workers=2,
        cluster_name='run vep',
    )
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
