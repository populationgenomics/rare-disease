#!/usr/bin/env python3

"""
wraps the validation script(s)
"""
import os
from pathlib import Path
from argparse import ArgumentParser
from cloudpathlib import AnyPath

import hailtop.batch as hb
from cpg_utils.config import get_config

from cpg_utils.git import (
    prepare_git_job,
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
)
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    output_path,
    remote_tmpdir,
)


DEFAULT_IMAGE = get_config()['workflow']['driver_image']
assert DEFAULT_IMAGE

MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')
OUTPUT_VCF = output_path('variants_from_mt.vcf.gz')


def set_job_resources(
    job: hb.batch.job.Job,
    auth=False,
    git=False,
    image: str | None = None,
    prior_job: hb.batch.job.Job | None = None,
    memory: str = 'standard',
):
    """
    applied resources to the job
    :param job:
    :param auth: if true, authenticate gcloud in this container
    :param git: if true, pull this repository into container
    :param image:
    :param prior_job:
    :param memory:
    """
    # apply all settings
    job.cpu(2).image(image or DEFAULT_IMAGE).memory(memory).storage('20G')

    if prior_job is not None:
        job.depends_on(prior_job)

    if auth:
        authenticate_cloud_credentials_in_job(job)

    if git:
        # copy the relevant scripts into a Driver container instance
        prepare_git_job(
            job=job,
            organisation=get_organisation_name_from_current_directory(),
            repo_name=get_repo_name_from_current_directory(),
            commit=get_git_commit_ref_of_current_repository(),
        )


def mt_to_vcf(batch: hb.Batch, input_file: str, output_file: str):
    """
    takes a MT and converts to VCF
    :param batch:
    :param input_file:
    :param output_file:
    :return:
    """
    mt_to_vcf_job = batch.new_job(name='Convert MT to VCF')
    set_job_resources(mt_to_vcf_job, git=True, auth=True)

    job_cmd = (
        f'PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} '
        f'--input {input_file} '
        f'--output {output_file}'
    )
    copy_common_env(mt_to_vcf_job)
    mt_to_vcf_job.command(job_cmd)
    return mt_to_vcf_job


def main(input_file: str):
    """

    """

    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='run reanalysis (AIP)',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory='highmem',
    )

    input_path = Path(input_file)
    if input_path.suffix == '.mt':
        if AnyPath(OUTPUT_VCF).exists():
            print("no need to convert")
        else:
            # requires conversion to vcf
            mt_to_vcf(batch=batch, input_file=input_file, output_file=OUTPUT_VCF)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='input_path')
    args = parser.parse_args()
    main(input_file=args.i)

