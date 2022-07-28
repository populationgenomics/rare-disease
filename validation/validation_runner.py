#!/usr/bin/env python3

"""
wraps the validation script(s)
"""

import os
import sys
from pathlib import Path
from argparse import ArgumentParser

import hailtop.batch.job
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
QUERY_COMPARISON = os.path.join(os.path.dirname(__file__), "hail_query_validate.py")
OUTPUT_VCFS = output_path("single_sample_vcfs")
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


def mt_to_vcf(batch: hb.Batch, input_file: str, header_lines: str | None):
    """
    takes a MT and converts to VCF
    adds in extra header lines for VQSR filters
    :param batch:
    :param input_file:
    :param header_lines:
    :return:
    """
    mt_to_vcf_job = batch.new_job(name="Convert MT to VCF")

    copy_common_env(mt_to_vcf_job)

    set_job_resources(mt_to_vcf_job, git=True, auth=True)

    job_cmd = (
        f"PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} "
        f"--input {input_file} "
        f"--output {OUTPUT_VCFS} "
    )

    if header_lines:
        job_cmd += f"--additional_header {header_lines}"

    logger.info(f"Command MT>VCF: {job_cmd}")

    mt_to_vcf_job.command(job_cmd)
    return mt_to_vcf_job


def comparison_job(
    batch, sample, truth_vcf, bed: str, prior_job, annotate=False, combine=False
):
    """

    Parameters
    ----------
    batch :
    sample :
    truth_vcf :
    bed :
    prior_job :
    annotate :
    combine :

    Returns
    -------

    """
    if annotate and combine:
        print("please pick a single setting [annotate/combine/split(default)]")
        sys.exit(1)

    vcf = os.path.join(OUTPUT_VCFS, f"{sample}.vcf.bgz")

    job = batch.new_job(name=f"compare_{sample}")
    job.image(HAPPY_IMAGE)
    job.memory("20Gi")
    vcf_input = batch.read_input_group(**{"vcf": vcf, "index": vcf + ".tbi"})
    truth_input = batch.read_input_group(
        **{"vcf": truth_vcf, "index": truth_vcf + ".tbi"}
    )
    truth_bed = batch.read_input(bed)
    refgenome = batch.read_input(
        "gs://cpg-reference/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta"
    )

    if annotate:
        job.declare_resource_group(
            output={
                "calls.vcf.gz": "{root}/calls.vcf.gz",
                "calls.vcf.gz.tbi": "{root}/calls.vcf.gz.tbi",
                "baseline.vcf.gz": "{root}/baseline.vcf.gz",
                "baseline.vcf.gz.tbi": "{root}/baseline.vcf.gz.tbi",
                "summary.txt": "{root}/summary.txt",
            }
        )

    elif combine:
        job.declare_resource_group(
            output={
                "combined_output.vcf.gz": "{root}/output.vcf.gz",
                "combined_output.vcf.gz.tbi": "{root}/output.vcf.gz.tbi",
                "summary.txt": "{root}/summary.txt",
            }
        )

    else:
        job.declare_resource_group(
            output={
                "fp.vcf.gz": "{root}/fp.vcf.gz",
                "fp.vcf.gz.tbi": "{root}/fp.vcf.gz.tbi",
                "fn.vcf.gz": "{root}/fn.vcf.gz",
                "fn.vcf.gz.tbi": "{root}/fn.vcf.gz.tbi",
                "tp.vcf.gz": "{root}/tp.vcf.gz",
                "tp.vcf.gz.tbi": "{root}/tp.vcf.gz.tbi",
                "summary.txt": "{root}/summary.txt",
            }
        )

    # default mode is split
    mode = ""
    if annotate:
        mode = "-m annotate "
    if combine:
        mode = "-m combine "

    # in future don't regenerate SDF
    job_cmd = (
        f"java -jar -Xmx16G /vcfeval/RTG.jar format -o refgenome_sdf {refgenome} && "
        f"mv {vcf_input['vcf']} input.vcf.gz && "
        f"mv {vcf_input['index']} input.vcf.gz.tbi && "
        f"java -jar -Xmx16G /vcfeval/RTG.jar vcfeval "
        f"--decompose "
        f"{mode}"
        f"-t refgenome_sdf "
        f"-b {truth_input['vcf']} "
        f"-c input.vcf.gz "
        f"-o {job.output} "
        f"--bed-regions={truth_bed} && "
        f"ls {job.output}"
    )

    job.command(job_cmd)
    job.depends_on(prior_job)
    batch.write_output(job.output, os.path.join(output_path("comparison"), sample))
    return job


def run_hail_benchmark(prior_job):
    pass


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
        # requires conversion to vcf
        prior_job = mt_to_vcf(
            batch=batch,
            input_file=input_file,
            header_lines=header,
        )

    # compare syndip
    _comparison_job = comparison_job(
        batch=batch,
        sample="syndip",
        prior_job=prior_job,
        bed="gs://cpg-reference/validation/syndip/regions/syndip.b38_20180222.bed",
        truth_vcf="gs://cpg-validation-test/syndip/syndip_truth.vcf.gz",
        combine=True,
    )
    # compare hg001
    _comparison_job = comparison_job(
        batch=batch,
        sample="na12878_kccg",
        prior_job=prior_job,
        bed="gs://cpg-validation-test/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
        truth_vcf="gs://cpg-reference/validation/giab/truth/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        combine=True,
    )

    batch.run(wait=False)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", help="input_path")
    parser.add_argument("--header", help="header_lines_file", default=None)
    args = parser.parse_args()
    main(input_file=args.i, header=args.header)
