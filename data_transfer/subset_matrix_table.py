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
import sys

from argparse import ArgumentParser

import hail as hl
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    output_path,
    query_command,
    remote_tmpdir,
)


def subset_to_samples(matrix: hl.MatrixTable, samples: list[str]) -> hl.MatrixTable:
    """
    reduce the MatrixTable to a subset
    Parameters
    ----------
    matrix :
    samples :

    Returns
    -------

    """

    mt_samples = matrix.s.collect()

    # check the subset exists
    missing_samples = [sam for sam in samples if sam not in mt_samples]
    if missing_samples:
        raise Exception(f"Sample(s) missing from subset: {','.join(missing_samples)}")

    filt_mt = matrix.filter_cols(matrix.s in samples)

    # # optional - filter to variants with at least one alt call in these samples
    # call_filt_matrix = filt_mt.filter_entries(filt_mt.GT.is_non_ref())

    return filt_mt


def subset_to_locus(matrix: hl.MatrixTable, chrom: str, pos: int) -> hl.MatrixTable:
    """

    Parameters
    ----------
    matrix :
    chrom :
    pos :

    Returns
    -------

    """

    locus = hl.Locus(contig=chrom, position=pos, reference_genome="GRCh38")
    matrix = matrix.filter_rows(matrix.locus == locus)
    if matrix.count_rows() == 0:
        raise Exception(f"No rows remain after applying Locus filter {locus}")
    return matrix


def main(
    mt_path: str,
    output_root: str,
    samples: list[str],
    vcf: bool,
    chrom: str | None,
    pos: int | None,
):
    """

    Parameters
    ----------
    mt_path : path to input MatrixTable
    output_root :
    samples :
    vcf :
    chrom :
    pos :

    Returns
    -------

    """

    hl.init(default_reference="GRCh38")
    matrix = hl.read_matrix_table(mt_path)

    if samples:
        matrix = subset_to_samples(matrix, samples)

    if chrom and pos:
        matrix = subset_to_locus(matrix=matrix, chrom=chrom, pos=pos)

    # write the MT to a new output path
    matrix.write(output_path(f"{output_root}.mt"))

    # if VCF, export as a VCF as well
    if vcf:
        vcf_output = output_path(f"{output_root}.vcf.bgz")
        hl.export_vcf(matrix, vcf_output, tabix=True)


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
    subset_job = batch.new_python_job("run matrix subsetting")
    subset_job.image(get_config()["workflow"]["driver_image"])
    # authenticate_cloud_credentials_in_job(subset_job)
    subset_job.call(
        main,
        args.i,
        args.out,
        args.s,
        args.vcf,
        args.chr,
        args.pos,
    )
    batch.run(wait=False)
