"""
Single purpose script to export a VCF version of a HT. This is currently specific to the [REDACTED] use case, but
could be made generic.

If this became generic, we would need to take a list/config for the fields which need to be migrated to the vcf INFO.
"""

import argparse

import hail as hl
from cpg_utils import hail_batch


def write_ht_as_vcf(ht_path: str, output_path: str) -> None:
    """Reads the HT which was just written, moves annotation into INFO, and writes as a VCF."""
    ht = hl.read_table(ht_path)
    ht = ht.transmute(
        info=hl.Struct(
            avis=ht.avis,
        ),
    )
    hl.export_vcf(ht, output_path, tabix=True)


if __name__ == '__main__':
    hail_batch.init_batch()

    parser = argparse.ArgumentParser(description='Read a VCF and export it as a VCF.')
    parser.add_argument(
        '--input',
        required=True,
        help='HT Path to write as VCF',
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output path for VCF of Hail Table',
    )
    args = parser.parse_args()

    write_ht_as_vcf(args.output, args.output)
