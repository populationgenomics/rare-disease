#!/usr/bin/env python3

"""
Takes a path to a MatrixTable and a name prefix
Writes data into the test bucket for the dataset

Optionally sample IDs and a locus can be provided to reduce the output
If sample IDs are specified, output subset contains only those
--format arg sets the output type (mt|vcf|both), default is MT
Optionally supply both --chr and --pos; subset to a specific locus
The pos format can be a single int, or a "start-end"
Default behaviour is to remove any sites where the requested samples
    are all HomRef. All sites can be retained with --keep_ref
"""

from argparse import ArgumentParser
import logging
import sys

import hail as hl

from cpg_utils.hail_batch import output_path, init_batch
from cpg_utils.config import get_config


def subset_to_samples(
    mt: hl.MatrixTable, samples: set[str], keep_hom_ref: bool
) -> hl.MatrixTable:
    """
    checks the requested sample subset exists in this joint call
    reduces the MatrixTable to a sample subset

    Parameters
    ----------
    mt : the current MatrixTable object
    samples : a set of samples requested for the subset selection
    keep_hom_ref : if False, remove sites where no remaining samples have alt calls

    Returns
    -------
    The original MT reduced to only the selected samples, and depending on the
    keep_hom_ref flag, only sites where that sample subset contain alt calls

    """

    missing_samples = samples.difference(set(mt.s.collect()))
    if missing_samples:
        raise Exception(f'Sample(s) missing from subset: {", ".join(missing_samples)}')

    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # optional - filter to variants with at least one alt call in these samples
    if not keep_hom_ref:
        mt = hl.variant_qc(mt)
        mt.filter_rows(mt.variant_qc.n_non_ref > 0)
        mt = mt.drop('variant_qc')

    return mt


def subset_to_locus(mt: hl.MatrixTable, locus: hl.IntervalExpression) -> hl.MatrixTable:
    """
    Subset the provided MT to a single locus - fail if the variant is absent

    Parameters
    ----------
    mt : the current MatrixTable object
    locus : a hail LocusExpression indicating the range of positions to select

    Returns
    -------
    The subset of the MatrixTable overlapping the indicated locus
    """

    mt = mt.filter_rows(locus.contains(mt.locus))
    if mt.count_rows() == 0:
        raise Exception(f'No rows remain after applying Locus filter {locus}')
    return mt


def main(
    mt_path: str,
    output_root: str,
    samples: set[str],
    out_format: str,
    locus: hl.IntervalExpression | None,
    keep_hom_ref: bool,
):
    """

    Parameters
    ----------
    mt_path : path to input MatrixTable
    output_root : prefix for file naming
    samples : a set of samples to reduce the joint-call to
    out_format : whether to write as a MT, VCF, or Both
    locus : an optional parsed interval for locus-based selection
    keep_hom_ref : if true, retain all sites in the subset

    Returns
    -------

    """

    mt = hl.read_matrix_table(mt_path)

    if samples:
        mt = subset_to_samples(mt, samples=samples, keep_hom_ref=keep_hom_ref)

    if isinstance(locus, hl.IntervalExpression):
        mt = subset_to_locus(mt, locus=locus)

    # create the output path; make sure we're only ever writing to test
    actual_output_path = output_path(output_root).replace(
        f'cpg-{get_config()["workflow"]["dataset"]}-main',
        f'cpg-{get_config()["workflow"]["dataset"]}-test',
    )

    if out_format in ['mt', 'both']:
        # write the MT to a new output path
        mt.write(f'{actual_output_path}.mt', overwrite=True)

    # if VCF, export as a VCF as well
    if out_format in ['vcf', 'both']:
        hl.export_vcf(mt, f'{actual_output_path}.vcf.bgz', tabix=True)


def clean_locus(contig: str, pos: str) -> hl.IntervalExpression | None:
    """

    Parameters
    ----------
    contig : the contig string, e.g. 'chr4'
    pos : either a single coordinate or two in the form start-stop

    Returns
    -------
    A parsed hail locus. For a point change this will be the result of
    parsing "contig:pos-pos+1"
    """
    if not all([contig, pos]):
        return None

    if '-' in pos:
        start, end = map(int, pos.split('-'))

        # quick validation that we only received 2 values
        assert isinstance(start, int)
        assert isinstance(end, int)

        assert start <= end

        if start == end:
            end += 1
    else:
        start = int(pos)
        end = start + 1

    return hl.parse_locus_interval(f'{contig}:{start}-{end}')


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the input MatrixTable', required=True)
    parser.add_argument(
        '--out',
        help='Prefix for MT/VCF name\n'
        '("output" will become output.vcf.bgz or output.mt)',
        required=True,
    )
    parser.add_argument(
        '-s', help='One or more sample IDs, whitespace delimited', nargs='+'
    )
    parser.add_argument(
        '--format',
        help='Write output in this format',
        default='mt',
        choices=['both', 'mt', 'vcf'],
    )
    parser.add_argument('--chr', help='Contig portion of a locus', required=False)
    parser.add_argument(
        '--pos',
        help='Pos portion of a locus. Can be "12345" or "12345-67890" for a range',
        required=False,
    )
    parser.add_argument(
        '--keep_ref',
        help='Output will retain all sites, even where the sample subset is HomRef',
        action='store_true',
    )
    args, unknown = parser.parse_known_args()

    if unknown:
        raise Exception(f'Unknown args, could not parse: "{unknown}"')

    if any([args.chr, args.pos]) and not all([args.chr, args.pos]):
        raise Exception(
            f'When defining a Locus, provide both Chr & Pos: {args.chr}, {args.pos}'
        )

    init_batch()
    locus_interval = clean_locus(args.chr, args.pos)

    main(
        mt_path=args.i,
        output_root=args.out,
        samples=set(args.s),
        out_format=args.format,
        locus=locus_interval,
        keep_hom_ref=args.keep_ref,
    )
