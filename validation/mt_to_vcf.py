"""
Takes an input MT, and extracts a VCF-format representation.
"""
import os
from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch
from cloudpathlib import AnyPath


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

    # extract all present samples into a separate file
    for sample in matrix.s.collect():

        sample_path = os.path.join(output_path, f"{sample}.vcf.bgz")

        if AnyPath(sample_path).exists():
            print(f"no action taken, {sample_path} already exists")
            continue

        # filter to this column
        ss_matrix = matrix.filter_cols(matrix.s == sample)

        # filter to this sample's genotype calls
        ss_matrix = ss_matrix.filter_entries(ss_matrix.GT.is_non_ref())

        hl.export_vcf(
            ss_matrix,
            sample_path,
            append_to_header=additional_header,
            tabix=True,
        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        help="input MatrixTable path",
    )
    parser.add_argument("--output", type=str, help="directory to write VCFs out to")
    parser.add_argument(
        "--additional_header",
        type=str,
        help="path to file containing any additional header lines",
        default=None,
    )
    args = parser.parse_args()
    main(
        input_mt=args.input,
        output_path=args.output,
        additional_header=args.additional_header,
    )
