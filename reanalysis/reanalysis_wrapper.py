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
from analysis_runner.git import (
    prepare_git_job,
    get_repo_name_from_current_directory,
    get_git_commit_ref_of_current_repository,
)

import click


# used to provide all VEP105 consequences, silences Slivar errors
PANELAPP_JSON_OUT = output_path('panelapp_137_data.json')
HAIL_VCF_OUT = output_path("hail_classified.vcf.bgz")
COMP_HET_VCF_OUT = output_path("hail_comp_het.vcf.bgz")

# location of the Slivar Docker image
AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
SLIVAR_TAG = 'slivar:v0.2.7'
SLIVAR_IMAGE = f'{AR_REPO}/{SLIVAR_TAG}'


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
@click.option('--ped', 'ped_file', help='ped file for this analysis')
def main(matrix_path: str, panelapp_date: str, config_json: str, ped_file: str):
    """
    Description
    :param matrix_path:
    :param panelapp_date:
    :param config_json:
    :param ped_file:
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(
        name='run_reanalysis', backend=service_backend, cancel_after_n_failures=1
    )

    # panelapp and hail script paths
    panelapp_script = os.path.join(os.path.dirname(__file__), 'panelapp_extraction.py')
    hail_script = os.path.join(os.path.dirname(__file__), 'hail_classifier.py')

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #
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

    # copy the relevant scripts into an instance
    # of the Driver container
    prepare_git_job(
        job=panelapp_job,
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )
    panelapp_job.image(os.getenv('DRIVER_IMAGE'))

    # retrieve the output file, writing to the output bucket
    batch.write_output(panelapp_job.panel_json, PANELAPP_JSON_OUT)

    # ----------------------- #
    # run hail classification #
    # ----------------------- #
    # insert skip clause here if output already exists
    hail_job = dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        script=f'python3 {hail_script} '
        f'--mt {matrix_path} '
        f'--pap {PANELAPP_JSON_OUT} '
        f'--config {config_json} '
        f'--output {HAIL_VCF_OUT}',
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
    # required?
    set_job_resources(hail_job)

    # ------------------------------------------- #
    # slivar compound het check & class 4 removal #
    # ------------------------------------------- #

    # set up the slivar job
    slivar_job = batch.new_job(name='slivar_reanalysis_stage')
    set_job_resources(slivar_job)
    slivar_job.image(SLIVAR_IMAGE)

    # copy in VCF from step 2 as input (HAIL_VCF_OUT)
    # copy in PED file for cohort
    # run tabix on the input resource file
    # run comp-het discovery on the file
    # this creates a new VCF, exclusively containing comp-hets
    slivar_job.command(
        (
            'tabix -p vcf input_vcf_resource; '
            'CSQ_FIELD="COMPOUND_CSQ" slivar compound-hets '
            f'--ped {batch.read_input(ped_file)} '
            f'-v {batch.read_input(HAIL_VCF_OUT)} | '
            f'bgzip -c -@ 4 > {slivar_job.out_vcf};'
        )
    )
    batch.write_output(slivar_job.out_vcf, COMP_HET_VCF_OUT)

    # run the batch
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
