from cpg_utils import hail_batch
import hail as hl
import argparse
# 1. Initialize Hail with the Service Backend for Batch
# This assumes you have your environment variables/billing project configured
hail_batch.init_batch()
parser = argparse.ArgumentParser(description='Import and merge Hail Tables')
parser.add_argument('--input_file1', type=str, help='Path to the input hail table')
parser.add_argument('--input_file2', type=str, help='Path to the input hail table')
parser.add_argument('--output', '-o', type=str, default=None,
                    help='Output path for merged Hail table (default: merged.ht)')

args = parser.parse_args()
MT_A = args.input_file1
MT_B = args.input_file2
DESTINATION = args.output if args.output else 'merged.mt'

def merge_matrix_tables(path_a, path_b, output_path):
    # Load the two Matrix Tables
    mt1 = hl.read_matrix_table(path_a)
    mt2 = hl.read_matrix_table(path_b)

    # 2. Perform the union
    merged_mt = mt1.union_cols(mt2)

    # 3. Write the result to a persistent location
    merged_mt.write(output_path, overwrite=True)

    print(f"Merge complete. Final MT saved to: {output_path}")


merge_matrix_tables(MT_A, MT_B, DESTINATION)
