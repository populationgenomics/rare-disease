#!/usr/bin/env python3

"""
script to run the whole multi-file AVI data ingestion in one process
requires that all component files have the required header (see reheader_files.sh)

designed to operate in QOB, with no localisation of TSV/tables

python3 .../full_avi_ingestion.py \
    /path/to/tsv1 /path/to/tsv2 /path/to/tsv3 \
    -o /path/to/output.ht \
    -t /path/to/tmp/bucket
"""

import argparse

import hail as hl
from cpg_utils import Path, hail_batch, to_path

REF_GENOME = "GRCh38"

# Define the input types for the initial table import
# We import chrom/pos/ref/alt as strings/ints first to transform them
INPUT_TYPES = {
    "#CHROM": hl.tstr,
    "POS": hl.tint32,
    "REF": hl.tstr,
    "ALT": hl.tstr,
    "avis": hl.tfloat64,
    "phred": hl.tfloat64,
}


def write_single_table(table_path: str, tmp_dir: Path) -> hl.Table:
    """
    Runs a transformation and checkpoint of a TSV into a Hail table

    Returns:
        the checkpointed, transformed, table
    """

    # this is a location in temp, not the final location. Create a name based on this TSV's filename
    ht_out = tmp_dir / f"{to_path(table_path).stem}.ht"

    # don't repeat work if we already made this table
    if (ht_out / "_SUCCESS").exists():
        return hl.read_table(str(ht_out))

    # load the TSV
    ht = hl.import_table(table_path, types=INPUT_TYPES, delimiter="\t", force_bgz=True)

    # Apply sigmoid transformation to avis column
    ht = ht.annotate(normalised_avis=hl.expit(ht.avis))
    ht = ht.rename({"avis": "raw_avis"})

    # create locus & alleles to key on, matching the compound key in our MatrixTables
    ht = ht.transmute(
        locus=hl.locus(ht["#CHROM"], ht.POS, reference_genome=REF_GENOME),
        alleles=[ht.REF, ht.ALT],
    )

    # 4. Key the table by locus and alleles
    ht = ht.key_by("locus", "alleles")

    # todo extend this with quantile when it's available across all component tables
    ht = ht.select("raw_avis", "normalised_avis", "phred")

    # Write the table to disk in Hail format for later use
    return ht.checkpoint(str(ht_out))


def merge_hail_tables(hail_tables: list[hl.Table], output_path: str) -> None:
    """Merge an arbitrary number of Hail Tables, assumes matching schema and non-overlapping rows."""

    # perform the multi-table union (concatenate rows)
    merged_ht = hail_tables[0].union(*hail_tables[1:])

    # write result to a persistent location
    merged_ht.write(output_path, overwrite=True)

    print(f"Merge complete. Final HT saved to: {output_path}")


def write_ht_as_vcf(ht_path: str, output_path: str) -> None:
    """Reads the HT which was just written, moves annotation into INFO, and writes as a VCF."""
    ht = hl.read_table(ht_path)
    ht = ht.transmute(
        info=hl.Struct(
            raw_avis=ht.raw_avis,
            normalised_avis=ht.normalised_avis,
            phred_avis=ht.phred,
        ),
    )
    hl.export_vcf(ht, output_path, tabix=True)

    print(f"Extraction complete, VCF representation saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tsvs", nargs="+")
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("-t", "--tmp", required=True)
    args = parser.parse_args()

    # start a QOB instance
    hail_batch.init_batch()

    assert args.output_path.endswith(".ht")

    tmp_as_path = to_path(args.tmp)

    # run transformations on each table
    all_tables = [
        write_single_table(table_path=table_path, tmp_dir=tmp_as_path)
        for table_path in args.tsvs
    ]

    # merge all tables into a single table
    merge_hail_tables(all_tables, args.output_path)

    # also write as a VCF (for echtvar transformation)
    write_ht_as_vcf(args.output_path, args.output_path.replace(".ht", ".vcf.bgz"))
