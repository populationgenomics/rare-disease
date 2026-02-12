import argparse

import hail as hl
from cpg_utils import hail_batch

# 1. Initialize Hail with the Service Backend for Batch
hail_batch.init_batch()

parser = argparse.ArgumentParser(description='Import and merge Hail Tables')
parser.add_argument(
    '--input',
    required=True,
    nargs='+',
    help='Paths to the input Hail tables to merge',
)
parser.add_argument(
    '--output',
    required=True,
    help='Output path for merged Hail Table',
)
parser.add_argument(
    '--vcf',
    required=False,
    default=None,
    help='Optional output path for the same data as a VCF',
)

args = parser.parse_args()


def merge_hail_tables(table_paths: list[str], output_path: str) -> None:
    """Merge an arbitrary number of Hail Tables, assumes matching schema and non-overlapping rows."""

    assert len(table_paths) > 1

    opened_tables = [hl.read_table(each_path) for each_path in table_paths]

    # 2. Perform the union (concatenate rows)
    merged_ht = opened_tables[0].union(*opened_tables[1:])

    # 3. Write the result to a persistent location
    merged_ht.write(output_path, overwrite=True)

    print(f"Merge complete. Final HT saved to: {output_path}")


def write_ht_as_vcf(ht_path: str, output_path: str) -> None:
    """Reads the HT which was just written, moves annotation into INFO, and writes as a VCF."""
    ht = hl.read_table(ht_path)
    ht = ht.transmute(
        info=hl.Struct(raw_avis=ht.raw_avis, normalised_avis=ht.normalised_avis),
    )
    hl.export_vcf(ht, output_path, tabix=True)


merge_hail_tables(table_paths=args.input, output_path=args.output)

# optional second write as a VCF
if args.vcf:
    write_ht_as_vcf(args.output, args.vcf)
