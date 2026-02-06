import argparse

import hail as hl
from cpg_utils import hail_batch

hail_batch.init_batch()

# Parse command line arguments
parser = argparse.ArgumentParser(description='Import TSV file into Hail Table format')
parser.add_argument(
    '--input',
    help='Path to the input TSV file (can be gzipped)',
    required=True,
)
parser.add_argument(
    '--output',
    help='Output path for Hail table',
    required=True,
)
args = parser.parse_args()

reference_genome = 'GRCh38'

# 1. Define the input types for the initial table import
# We import chrom/pos/ref/alt as strings/ints first to transform them
input_types = {
    '#CHROM': hl.tstr,
    'POS': hl.tint32,
    'REF': hl.tstr,
    'ALT': hl.tstr,
    'avis': hl.tfloat64,
}

# 2. Load the TSV
ht = hl.import_table(args.input, types=input_types, delimiter='\t', force_bgz=True)

# 3. Transform to standard Hail genomic format
# We create a 'locus' object and an 'alleles' array
ht = ht.transmute(
    locus=hl.locus(ht['#CHROM'], ht.POS, reference_genome=reference_genome),
    alleles=[ht.REF, ht.ALT],
)

# 4. Key the table by locus and alleles
ht = ht.key_by('locus', 'alleles')

ht.describe()

# Write the table to disk in Hail format for later use
ht.write(args.output, overwrite=True)
print(f'Table successfully written to: {args.output}')
