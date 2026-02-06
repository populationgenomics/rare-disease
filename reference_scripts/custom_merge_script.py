import argparse

import hail as hl
from cpg_utils import hail_batch

# 1. Initialize Hail with the Service Backend for Batch
# This assumes you have your environment variables/billing project configured
hail_batch.init_batch()
parser = argparse.ArgumentParser(description='Import and merge Hail Tables')
parser.add_argument(
    '--input_file1',
    type=str,
    required=True,
    help='Path to the input hail table',
)
parser.add_argument(
    '--input_file2',
    type=str,
    required=True,
    help='Path to the input hail table',
)
parser.add_argument(
    '--output',
    '-o',
    type=str,
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
HT_A = args.input_file1
HT_B = args.input_file2
DESTINATION = args.output


def merge_hail_tables(path_a: str, path_b: str, output_path: str) -> None:
    # Load the two Hail Tables
    ht1 = hl.read_table(path_a)
    ht2 = hl.read_table(path_b)

    # 2. Perform the union (concatenate rows)
    merged_ht = ht1.union(ht2)

    # 3. Write the result to a persistent location
    merged_ht.write(output_path, overwrite=True)

    print(f"Merge complete. Final HT saved to: {output_path}")


def write_ht_as_vcf(ht_path: str, output_path: str) -> None:
    """Reads the HT which was just written, moves annotation into INFO, and writes as a VCF."""
    ht = hl.read_table(ht_path)
    ht = ht.transmute(
        info=hl.Struct(
            avis=ht.avis,
        ),
    )
    hl.export_vcf(ht, output_path, tabix=True)

merge_hail_tables(HT_A, HT_B, DESTINATION)

# optional second write as a VCF
if args.vcf:
    write_ht_as_vcf(DESTINATION, args.vcf)
