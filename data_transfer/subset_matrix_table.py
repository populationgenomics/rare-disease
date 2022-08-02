#!/usr/bin/env python3

"""
Takes a path to a MatrixTable, and optionally sample IDs
Reads in the MT from provided input path
If sample IDs are specified, subset containing only those
Write the [reduced] MT into Test
Optionally with the --vcf flag, write as a VCF
Optionally with both --chr and --pos specified, subset to a specific locus
"""

from argparse import ArgumentParser
import logging
import sys

import hail as hl

from cpg_utils.hail_batch import output_path, init_batch
from cpg_utils.config import get_config


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

    filt_mt = matrix.filter_cols(hl.literal(set(samples)).contains(matrix.s))

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
    init_batch()
    matrix = hl.read_matrix_table(mt_path)

    if samples:
        matrix = subset_to_samples(matrix, samples)

    if chrom and pos:
        matrix = subset_to_locus(matrix=matrix, chrom=chrom, pos=pos)

    # create the output path; make sure we're only ever writing to test
    actual_output_path = output_path(output_root).replace(
        f'cpg-{get_config()["workflow"]["dataset"]}-main',
        f'cpg-{get_config()["workflow"]["dataset"]}-test',
    )

    # write the MT to a new output path
    matrix.write(f"{actual_output_path}.mt")

    # if VCF, export as a VCF as well
    if vcf:
        hl.export_vcf(matrix, f"{actual_output_path}.vcf.bgz", tabix=True)


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

    main(
        mt_path=args.i,
        output_root=args.out,
        samples=args.s,
        vcf=args.vcf,
        chrom=args.chr,
        pos=args.pos,
    )