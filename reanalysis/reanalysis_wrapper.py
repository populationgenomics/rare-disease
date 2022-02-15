#!/usr/bin/env python3


"""
wrapper for reanalysis process
"""


from typing import Union
import logging
import os
import hailtop.batch as hb
import hailtop.batch.job
from analysis_runner import dataproc, output_path

import click


def set_job_resources(job: Union[hailtop.batch.job.BashJob, hailtop.batch.job.Job]):
    """
    returns a new, appropriately resourced job
    :param job:
    """
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')


@click.command()
@click.option('--matrix', 'matrix_path', help='variant matrixtable to analyse')
@click.option('--pap_date', 'panelapp_date', help='panelapp date threshold')
@click.option(
    '--conf',
    'config_json',
    help='dictionary of runtime thresholds',
    required=False,
    default='gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json',
)
def main(matrix_path: str, panelapp_date: str, config_json: str):
    """
    Description
    :param matrix_path:
    :param panelapp_date:
    :param config_json:
    :return:
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(name='run_reanalysis', backend=service_backend)

    # panelapp and hail script paths
    panelapp_script = os.path.join(os.path.dirname(__file__), 'panelapp_extraction.py')
    hail_script = os.path.join(os.path.dirname(__file__), 'hail_classifier.py')

    # grab panelapp content in JSON form
    panelapp_job = batch.new_job(name='parse panelapp')
    set_job_resources(panelapp_job)
    panelapp_command = (
        f'python3 {panelapp_script} '
        f'--id 137 '
        f'--out {panelapp_job.panel_json} '
        f'--date {panelapp_date}'
    )
    logging.info('PanelApp Process trigger: %s', panelapp_command)
    panelapp_job.command(panelapp_command)
    batch.write_output(panelapp_job.panel_json, output_path('panelapp_137_data.json'))

    # don't run this for now
    if panelapp_date == '':
        # create a hail job
        hail_job = dataproc.hail_dataproc_job(
            batch=batch,
            worker_machine_type='n1-highmem-8',
            worker_boot_disk_size=200,
            secondary_worker_boot_disk_size=200,
            script=f'{hail_script} '
            f'--mt {matrix_path} '
            f'--pap {output_path("panelapp_137_data.json")} '
            f'--config {config_json} '
            f'--output {output_path("hail_classified.vcf.bgz")}',
            max_age='8h',
            init=[
                'gs://cpg-reference/hail_dataproc/install_common.sh',
                'gs://cpg-reference/vep/vep-GRCh38.sh',  # install and configure VEP 105
            ],
            job_name='run_vep',
            num_secondary_workers=10,
            num_workers=2,
            cluster_name='run vep',
        )
        set_job_resources(hail_job)

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
