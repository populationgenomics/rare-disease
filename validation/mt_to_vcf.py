"""
Takes an input MT, and extracts a VCF-format representation.
"""


import logging
import os
from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch
from cloudpathlib import AnyPath


def main(
    input_mt: str,
    validation_samples: set[str],
    output_path: str,
    additional_header: str | None = None,
):
    """
    takes an input MT, and reads it out as a VCF
    :param input_mt:
    :param validation_samples: all samples to target
    :param output_path:
    :param additional_header: file containing lines to append to header
    :return:
    """

    init_batch()

    full_mt = hl.read_matrix_table(input_mt)

    samples_in_jc = set(full_mt.s.collect()).intersection(validation_samples)
    logging.info(
        f'Extracting these samples from the joint-call: {" ".join(samples_in_jc)}'
    )

    # extract all present samples into a separate file
    for sample in samples_in_jc:

        logging.info(f'Processing: {sample}')

        sample_path = os.path.join(output_path, f'{sample}.vcf.bgz')

        if AnyPath(sample_path).exists():
            print(f'no action taken, {sample_path} already exists')
            continue

        # filter to this column
        mt = full_mt.filter_cols(full_mt.s == sample)

        # filter to this sample's non-ref calls
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)
        hl.export_vcf(
            mt,
            sample_path,
            append_to_header=additional_header,
            tabix=True,
        )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        type=str,
        help='input MatrixTable path',
    )
    parser.add_argument('--output', type=str, help='directory to write VCFs out to')
    parser.add_argument(
        '--additional_header',
        type=str,
        help='path to file containing any additional header lines',
        default=None,
    )
    parser.add_argument('--samples', nargs='+', help='validation samples to target')
    args = parser.parse_args()
    main(
        input_mt=args.input,
        output_path=args.output,
        additional_header=args.additional_header,
        validation_samples=set(args.samples),
    )
