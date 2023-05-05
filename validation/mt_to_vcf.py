"""
Takes an input MT, and extracts a VCF-format representation.
"""


from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch


def main(input_mt: str, sample: str, write_path: str):
    """
    takes an input MT, and reads it out as a VCF
    :param input_mt:
    :param sample: the sample to extract
    :param write_path:
    :return:
    """

    init_batch()

    # read full MT
    mt = hl.read_matrix_table(input_mt)

    # filter to this column
    mt = mt.filter_cols(mt.s == sample)

    # filter to this sample's non-ref calls
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    # filter out any Filter-failures
    mt = mt.filter_rows(mt.filters.length() == 0)

    hl.export_vcf(mt, write_path, tabix=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='input MT path')
    parser.add_argument('-o', help='directory to write VCFs out to')
    parser.add_argument('-s', help='validation sample to target')
    args = parser.parse_args()
    main(input_mt=args.i, write_path=args.o, sample=args.s)
