#!/usr/bin/env python3


"""
wrapper for reanalysis process
"""


from typing import Any, Dict, Union
import json
import logging
import os

from google.cloud import storage
import hailtop.batch as hb
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
REHEADERED_OUT = output_path("hail_classes_reheadered.vcf.bgz")
COMP_HET_VCF_OUT = output_path("hail_comp_het.vcf.bgz")
MT_OUT_PATH = output_path('hail_105_ac.mt')
RESULTS_HTML = output_path('summary_output.html')


# location of the Slivar Docker image
AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
SLIVAR_TAG = 'slivar:v0.2.7'
BCFTOOLS_TAG = 'bcftools:1.10.2--h4f4756c_2'
SLIVAR_IMAGE = f'{AR_REPO}/{SLIVAR_TAG}'
BCFTOOLS_IMAGE = f'{AR_REPO}/{BCFTOOLS_TAG}'

HAIL_SCRIPT = os.path.join(os.path.dirname(__file__), "hail_array_classifier.py")
PANELAPP_SCRIPT = os.path.join(os.path.dirname(__file__), "panelapp_extraction.py")
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), "validate_classifications.py")


def get_gcp_blob(bucket_path: str) -> storage.blob:
    """
    take a GCP bucket path to a file, read into a blob object
    :param bucket_path:
    :return: a blob representing the data
    """

    # split the full path to get the bucket and file path
    bucket = bucket_path.replace('gs://', '').split('/')[0]
    path = bucket_path.replace('gs://', '').split('/', maxsplit=1)[1]

    # create a client
    g_client = storage.Client()

    # obtain the blob of the data
    return g_client.get_bucket(bucket).get_blob(path)


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a GCP bucket path to a JSON file, read into an object
    this loop can read config files, or data
    :param bucket_path:
    :return:
    """

    json_blob = get_gcp_blob(bucket_path)

    # download_as_bytes method isn't available; but this returns bytes?
    return json.loads(json_blob.download_as_string())


def set_html_blob_metadata(filepath: str):
    """
    takes a path to an HTML file, and sets appropriate metadata
    """
    blob = get_gcp_blob(filepath)
    if blob is None or not blob.exists():
        logging.error('File %s doesn\'t exist, cannot set metadata', filepath)
        return

    # if some metadata is already populated, add/update
    if isinstance(blob.metadata, dict):
        blob.metadata['Content-Type'] = 'html'

    # otherwise create the dictionary
    else:
        blob.metadata = {'Content-Type': 'html'}

    # set the main ACL attribute (determines browser access)
    blob.acl.blob.content_type = 'html'
    blob.patch()

    logging.info('HTML metadata attributes have been added to %s', filepath)


def check_blob_exists(filepath: str) -> bool:
    """
    used to allow for the skipping of long running process
    if data already exists

    - for novel analysis runs, might need a force parameter
    if output folder is the same as a prev. run
    :param filepath:
    :return:
    """
    blob = get_gcp_blob(filepath)

    if blob is None:
        return False
    return blob.exists()


def set_job_resources(job: Union[hb.batch.job.BashJob, hb.batch.job.Job], git=False):
    """
    applied resources to the job
    :param job:
    :param git:
    """
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')
    if git:
        # copy the relevant scripts into a Driver container instance
        prepare_git_job(
            job=job,
            repo_name=get_repo_name_from_current_directory(),
            commit=get_git_commit_ref_of_current_repository(),
        )


def handle_panelapp_job(batch: hb.Batch, date: str) -> hb.batch.job.Job:
    """

    :param batch:
    :param date:
    :return:
    """
    panelapp_job = batch.new_job(name='parse panelapp')
    panelapp_job.image(os.getenv('DRIVER_IMAGE'))
    set_job_resources(panelapp_job, git=True)
    panelapp_command = (
        f'python3 {PANELAPP_SCRIPT} '
        f'--out {panelapp_job.panel_json} '
        f'--date {date}'
    )
    logging.info('PanelApp Process trigger: %s', panelapp_command)

    panelapp_job.command(panelapp_command)

    # retrieve the output file, writing to the output bucket
    batch.write_output(panelapp_job.panel_json, PANELAPP_JSON_OUT)
    return panelapp_job


def handle_hail_job(
    batch: hb.Batch, matrix: str, config: str, prior_job: hb.batch.job
) -> hb.batch.job.Job:
    """
    sets up the hail dataproc stage
    :param batch:
    :param matrix:
    :param config:
    :param prior_job:
    :return:
    """
    hail_job = dataproc.hail_dataproc_job(
        batch=batch,
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        script=f'{HAIL_SCRIPT} '
        f'--mt {matrix} '
        f'--pap {PANELAPP_JSON_OUT} '
        f'--config {config} '
        f'--output {HAIL_VCF_OUT} '
        f'--mt_out {MT_OUT_PATH}',
        max_age='8h',
        init=[
            'gs://cpg-reference/hail_dataproc/install_common.sh',
            'gs://cpg-reference/vep/vep-GRCh38.sh',  # install and configure VEP 105
        ],
        job_name='run hail reanalysis stage',
        num_secondary_workers=10,
        num_workers=2,
        cluster_name='hail_reanalysis_stage',
    )
    set_job_resources(hail_job)

    # don't start unless prior step is successful
    hail_job.depends_on(prior_job)

    return hail_job


def handle_reheader_job(
    batch: hb.Batch,
    local_vcf: str,
    prior_job: hb.batch.job,
    config_dict: Dict[str, Any],
) -> hb.batch.job:
    """
    runs the bcftools re-header process
    :param batch:
    :param local_vcf:
    :param prior_job:
    :param config_dict:
    :return:
    """

    bcftools_job = batch.new_job(name='bcftools_reheader_stage')
    set_job_resources(bcftools_job)
    bcftools_job.image(BCFTOOLS_IMAGE)

    bcftools_job.depends_on(prior_job)

    bcftools_job.declare_resource_group(
        vcf={'vcf': '{root}.vcf.bgz', 'vcf.tbi': '{root}.vcf.bgz.tbi'}
    )

    # reheader the VCF using BCFtools and sed
    desc = '##INFO=<ID=CSQ,Number=.,Type=String,Description="'

    conf_csq = config_dict.get('csq_string').replace('|', r'\|')
    new_format = rf"Format: '{conf_csq}'"

    bcftools_job.command(
        'set -ex; '
        f'bcftools view -h {local_vcf} | sed \'s/'
        f'{desc}">/{desc}{new_format}">/\' > new_header; '
        f'bcftools reheader -h new_header --threads 4 -o {bcftools_job.vcf["vcf"]} {local_vcf}; '
        f'tabix {bcftools_job.vcf["vcf"]}; '
    )
    return bcftools_job


def handle_slivar_job(
    batch: hb.Batch,
    reheadered_vcf: str,
    local_ped: str,
    prior_job: hb.batch.job,
) -> hb.batch.job:
    """
    set up the slivar job
    use PED and VCF as hail resource files
    tabix input resource file
    run comp-het discovery on the file
    creates new VCF  exclusively containing comp-hets

    Maybe do reheadering as a separate job? CSQ field could be broadly useful

    :param batch:
    :param reheadered_vcf:
    :param local_ped:
    :param prior_job:
    :return:
    """
    slivar_job = batch.new_job(name='slivar_reanalysis_stage')
    set_job_resources(slivar_job)
    slivar_job.image(SLIVAR_IMAGE)

    slivar_job.depends_on(prior_job)

    slivar_job.command(
        'export SLIVAR_QUIET="true"; '
        'slivar compound-hets '
        '--allow-non-trios '
        f'--ped {local_ped} '
        f'-v {reheadered_vcf} | '
        f'bgzip -c -@ 4 > {slivar_job.out_vcf};'
    )
    return slivar_job


@click.command()
@click.option('--matrix', 'matrix_path', help='variant matrix table to analyse')
@click.option('--pap_date', 'panelapp_date', help='panelapp date threshold')
@click.option(
    '--conf',
    'config_json',
    help='dictionary of runtime thresholds',
    required=False,
    default='gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json',
)
@click.option('--ped', 'ped_file', help='ped file for this analysis')
def main(
    matrix_path: str, panelapp_date: str, config_json: str, ped_file: str
):  # pylint: disable=too-many-locals

    """
    Description
    :param matrix_path:
    :param panelapp_date:
    :param config_json:
    :param ped_file:
    """

    config_dict = read_json_dict_from_path(config_json)

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(
        name='run_reanalysis', backend=service_backend, cancel_after_n_failures=1
    )

    # read ped and config files as a local batch resource
    # fiddle with the ped file so that we are _really_ doing singletons
    # SampleMetadata API to generate PED from *project*
    ped_in_batch = batch.read_input(ped_file)
    conf_in_batch = batch.read_input(config_json)

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #
    panelapp_job = handle_panelapp_job(batch=batch, date=panelapp_date)

    # ----------------------- #
    # run hail classification #
    # ----------------------- #
    # retain dependency flow if we skip the hail job
    prior_job = panelapp_job

    # complete absence of the VCF will cause an exception to be thrown
    # handled in method for now
    if not check_blob_exists(HAIL_VCF_OUT):
        hail_job = handle_hail_job(
            batch=batch,
            matrix=matrix_path,
            config=config_json,
            prior_job=panelapp_job,
        )
        prior_job = hail_job

    # copy the Hail output file into the remaining batch jobs
    hail_output_in_batch = batch.read_input_group(
        **{'vcf': HAIL_VCF_OUT, 'vcf.tbi': HAIL_VCF_OUT + '.tbi'}
    )

    # --------------------------------- #
    # bcftools re-headering of hail VCF #
    # --------------------------------- #
    bcftools_job = handle_reheader_job(
        batch=batch,
        local_vcf=hail_output_in_batch['vcf'],
        prior_job=prior_job,
        config_dict=config_dict,
    )

    batch.write_output(bcftools_job.vcf, REHEADERED_OUT)

    # ------------------------- #
    # slivar compound het check #
    # ------------------------- #
    slivar_job = handle_slivar_job(
        batch=batch,
        reheadered_vcf=bcftools_job.vcf['vcf'],
        local_ped=ped_in_batch,
        prior_job=prior_job,
    )

    batch.write_output(slivar_job.out_vcf, COMP_HET_VCF_OUT)

    # not yet in use??
    results_job = batch.new_job(name='finalise_results')
    # don't start unless prior jobs are successful
    results_job.depends_on(slivar_job)

    set_job_resources(results_job, git=True)

    # needs a container with either cyvcf2 or pyvcf inside
    # we could be gross here, and tuck in an installation?
    results_command = (
        'export MAMBA_ROOT_PREFIX="/root/micromamba" && '
        'micromamba install -y cyvcf2 --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge && '
        f'PYTHONPATH=$(pwd) python3 {RESULTS_SCRIPT} '
        f'--conf {conf_in_batch} '
        f'--class_vcf {bcftools_job.vcf["vcf"]} '
        f'--comp_het {slivar_job.out_vcf} '
        f'--pap {panelapp_job.panel_json} '
        f'--ped {ped_in_batch} '
        f'--out_path {results_job.ofile} '
    )
    logging.info('Results trigger: %s', results_command)

    results_job.command(results_command)
    results_job.image(os.getenv('DRIVER_IMAGE'))

    # write the final results file to an HTML
    batch.write_output(results_job.ofile, RESULTS_HTML)

    # # add a new job here just to overwrite the HTML metadata
    # metadata_job = batch.new_python_job()
    # metadata_job.call(set_html_blob_metadata, RESULTS_HTML)

    # run the batch
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
