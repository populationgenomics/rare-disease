"""
Script to extract a VCF representation of a Hail Table
input and output arguments are self-explanatory.
--fields takes a list of Strings, these represent fields in the HT schema which will be migrated to vcf.INFO.

This script doesn't accommodate nested fields.
"""

import argparse

import hail as hl
from cpg_utils import hail_batch


def write_ht_as_vcf(ht_path: str, output_path: str, fields: list[str]) -> None:
    """Reads the HT which was just written, moves annotation into INFO, and writes as a VCF."""
    ht = hl.read_table(ht_path)
    ht = ht.transmute(
        info=hl.Struct(**{f: ht[f] for f in fields}),
    )
    hl.export_vcf(ht, output_path, tabix=True)


if __name__ == '__main__':
    hail_batch.init_batch()

    parser = argparse.ArgumentParser(description='Read a HT and export it as a VCF.')
    parser.add_argument(
        '--input',
        required=True,
        help='Input Hail Table path to export as VCF',
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output path for VCF of Hail Table',
    )
    parser.add_argument(
        '--fields',
        nargs='+',
        help='Specific fields to move from the HT base schema to vcf.INFO.',
    )
    args = parser.parse_args()

    write_ht_as_vcf(args.input, args.output, fields=args.fields)
