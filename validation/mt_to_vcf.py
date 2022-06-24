"""
Takes an input MT, and extracts a VCF-format representation.
"""

from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch


def main(input_mt: str, output_path: str, additional_header: str | None = None):
    """
    takes an input MT, and reads it out as a VCF
    :param input_mt:
    :param output_path:
    :param additional_header: file containing lines to append to header
    :return:
    """
    init_batch()

    matrix = hl.read_matrix_table(input_mt)

    hl.export_vcf(
        matrix,
        output_path,
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
    parser.add_argument('--output', type=str, help='path to write VCF out to')
    parser.add_argument(
        '--additional_header',
        type=str,
        help='path to file containing any additional header lines',
        required=False,
        default=None,
    )
    args = parser.parse_args()
    main(
        input_mt=args.input,
        output_path=args.output,
        additional_header=args.additional_header,
    )
