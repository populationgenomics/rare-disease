#!/usr/bin/env python3

"""
wraps the validation script(s)
"""

import os
from pathlib import Path
from argparse import ArgumentParser

import hailtop.batch.job
from cloudpathlib import AnyPath
import logging

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
    image_path,
)


DEFAULT_IMAGE = get_config()["workflow"]["driver_image"]
HAPPY_IMAGE = image_path("happy-vcfeval")
assert DEFAULT_IMAGE, HAPPY_IMAGE

MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), "mt_to_vcf.py")
OUTPUT_VCF = output_path("variants_from_mt.vcf.bgz")

# create a logger
logger = logging.getLogger(__file__)


def set_job_resources(
    job: hb.batch.job.Job,
    auth=False,
    git=False,
    image: str | None = None,
    prior_job: hb.batch.job.Job | None = None,
    memory: str = "standard",
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
    job.cpu(2).image(image or DEFAULT_IMAGE).memory(memory).storage("20G")

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


def mt_to_vcf(
    batch: hb.Batch, input_file: str, output_file: str, header_lines: str | None
):
    """
    takes a MT and converts to VCF
    adds in extra header lines for VQSR filters
    :param batch:
    :param input_file:
    :param output_file:
    :param header_lines:
    :return:
    """
    mt_to_vcf_job = batch.new_job(name="Convert MT to VCF")

    copy_common_env(mt_to_vcf_job)

    set_job_resources(mt_to_vcf_job, git=True, auth=True)

    job_cmd = (
        f"PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} "
        f"--input {input_file} "
        f"--output {output_file} "
    )

    if header_lines:
        job_cmd += f"--additional_header {header_lines}"

    logger.info(f"Command MT>VCF: {job_cmd}")

    mt_to_vcf_job.command(job_cmd)
    return mt_to_vcf_job


def compare_syndip(
    batch: hailtop.batch.Batch, vcf: str, prior_job
) -> hailtop.batch.job.Job:
    """

    Parameters
    ----------
    batch : the batch to insert this job into
    vcf : input VCF file
    prior_job : a job to depend on, or None

    Returns
    -------

    """

    # this data should be supplied externally, but is hard coded for now
    syndip_truth = "gs://cpg-validation-test/syndip/syndip_truth.vcf.gz"
    syndip_bed = "gs://cpg-reference/validation/syndip/regions/syndip.b38_20180222.bed"
    syndip_sample = "SYNDIP"

    job = batch.new_job(name="compare_syndip")
    job.image(HAPPY_IMAGE)
    job.memory("20Gi")
    vcf_input = batch.read_input_group(**{"vcf": vcf, "index": vcf + ".tbi"})
    truth_input = batch.read_input_group(
        **{"vcf": syndip_truth, "index": syndip_truth + ".tbi"}
    )
    truth_bed = batch.read_input(syndip_bed)
    refgenome = batch.read_input(
        "gs://cpg-reference/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta"
    )

    job_cmd = (
        f"java -jar -Xmx16G /vcfeval/RTG.jar format -o refgenome_sdf {refgenome} && "
        f"mv {vcf_input['vcf']} > input.vcf.gz"
        f"mv {vcf_input['index']} > input.vcf.gz.tbi"
        f"java -jar -Xmx16G /vcfeval/RTG.jar vcfeval "
        f"-t refgenome_sdf "
        f"-b {truth_input['vcf']} "
        f"-c input.vcf.gz "
        f"-o {job.output} "
        f"--bed-regions={truth_bed} "
        f"--sample={syndip_sample} &&"
        f"ls {job.output}"
    )

    job.command(job_cmd)

    return prior_job


def main(input_file: str, header: str | None):
    """

    Parameters
    ----------
    input_file :
    header :

    Returns
    -------

    """

    service_backend = hb.ServiceBackend(
        billing_project=get_config()["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name="run validation bits and pieces",
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory="highmem",
    )

    prior_job = None

    input_path = Path(input_file)
    if input_path.suffix == ".mt":
        if AnyPath(OUTPUT_VCF).exists():
            print("no need to convert")
        else:
            # requires conversion to vcf
            prior_job = mt_to_vcf(
                batch=batch,
                input_file=input_file,
                output_file=OUTPUT_VCF,
                header_lines=header,
            )

    # now do some comparison-related things
    _prior_job = compare_syndip(batch=batch, vcf=OUTPUT_VCF, prior_job=prior_job)

    batch.run(wait=False)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", help="input_path")
    parser.add_argument("--header", help="header_lines_file", default=None)
    args = parser.parse_args()
    main(input_file=args.i, header=args.header)
