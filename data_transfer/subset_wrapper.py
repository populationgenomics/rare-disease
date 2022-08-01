#!/usr/bin/env python3

"""
Takes a path to a MatrixTable, and optionally sample IDs
Reads in the MT from provided input path
If sample IDs are specified, subset containing only those
Write the [reduced] MT into Test
Optionally with the --vcf flag, write as a VCF
Optionally with both --chr and --pos specified, subset to a specific locus
"""

import logging
import os
import sys

from argparse import ArgumentParser

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
    remote_tmpdir,
)

script_path = os.path.join(os.path.dirname(__file__), "subset_matrix_table.py")


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s",
        datefmt="%Y-%M-%d %H:%M:%S",
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument("-i", help="path to the input MatrixTable", required=True)
    parser.add_argument(
        "--out",
        help=(
            "prefix for the output MT or VCF (e.g. 'output' will become "
            "<output>/output.vcf.bgz or <output>/output.mt\n"
            "where <output> is determined by the analysis_runner -o setting"
        ),
        required=True,
    )
    parser.add_argument("-s", help="one or more sample IDs", nargs="+", default=[])
    parser.add_argument(
        "--vcf", help="write output as a VCF", action="store_true", default=False
    )
    parser.add_argument("--chr", help="chrom portion of a locus", required=False)
    parser.add_argument(
        "--pos", help="pos portion of a locus", required=False, type=int
    )
    args = parser.parse_args()

    if any([args.chr, args.pos]) and not all([args.chr, args.pos]):
        raise Exception(
            f"When defining a Locus, provide both Chr & Pos: {args.chr}, {args.pos}"
        )

    service_backend = hb.ServiceBackend(
        billing_project=get_config()["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name="AIP batch",
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory="highmem",
    )
    subset_job = batch.new_job("run matrix subsetting")
    subset_job.image(get_config()["workflow"]["driver_image"])
    authenticate_cloud_credentials_in_job(subset_job)
    prepare_git_job(
        job=subset_job,
        organisation=get_organisation_name_from_current_directory(),
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )

    sample_arg = f'-s {" ".join(args.s)}' if args.s else ""
    locus = f"--chr {args.chr} --pos {args.pos}" if args.chr else ""
    vcf_arg = "--vcf" if args.vcf else ""
    subset_job.command(
        f"python3 {script_path} -i {args.i} --out {args.out} {vcf_arg} {sample_arg} {locus}"
    )
    batch.run(wait=False)
