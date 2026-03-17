"""
Script to read a matrix table and reformat the annotations to be more discoverable

Goals:
 - read MT data in
 - reformat annotations to be more like a Table/TSV
 - remove unnecessary columns/fields
 - checkpoint the data after reduction

How:

 - migrate important nested fields up to the top level of the MT
    - retain the data, but in a more easily searchable way
 - aggregate all CPG IDs with variants up to the top level of the MT
    - allows us to remove all entry fields, whilst

Learn to love the Hail documentation, personally I think it's the best thing about Hail
https://hail.is/docs/0.2/expressions.html
"""

from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import init_batch


def load_in_mt(mt_path: str) -> hl.MatrixTable:
    """Read in the MT."""
    return hl.read_matrix_table(mt_path)


def main(input_path: str, output_path: str) -> None:
    """
    do all the things

    Args:
        input_path: where to find the input MT (product of annotate_cohort, MT)
        output_path: either a HT, or a TSV (see section below on writing output)
    """

    # start a Query-on-batch runtime
    init_batch(
        worker_cores=8,
        worker_memory='highmem',
        driver_cores=8,
        driver_memory='highmem',
    )

    # get the MatrixTable
    mt = load_in_mt(input_path)

    # 1. Count exactly how many rows are missing the 'avis' score
    missing_avis_count = mt.aggregate_rows(hl.agg.count_where(hl.is_missing(mt.avis)))
    print(f"Number of variants missing an AVI score: {missing_avis_count}")

    # 2. (Optional) Visually inspect what those missing rows look like
    print("Here is a peek at the rows with missing AVIS scores:")
    mt.filter_rows(hl.is_missing(mt.avis)).rows().select('avis').show(5)

    # keep a record of all the fields we want to keep in the final output (used in a dataset.select command to throw everything else away)
    fields_to_keep: list[str] = []

    # 1. shift some of the simple annotations up to the top level, and/or give them a simple name
    # some of this has already been done in the pipeline:
    # https://github.com/populationgenomics/cpg-flow-seqr-loader/blob/main/src/cpg_seqr_loader/scripts/annotate_cohort.py#L265-L292
    # but some of the fields remain nested. see schema in comments on https://cpg-populationanalysis.atlassian.net/browse/RD-757
    mt = mt.annotate_rows(
        primate_ai=mt.primate_ai.score,
        revel=mt.dbnsfp.REVEL_score,
        splice_ai_delta=mt.splice_ai.delta_score,
        cadd=mt.cadd.PHRED,
        eigen=mt.eigen.Eigen_phred,
        gnomad_genomes_faf=mt.gnomad_genomes.FAF_AF,
        gnomad_exomes_faf=mt.gnomad_exomes.FAF_AF,
    )

    # add these fields to the list to keep
    fields_to_keep.extend(
        [
            'primate_ai',
            'revel',
            'splice_ai_delta',
            'cadd',
            'eigen',
            'gnomad_genomes_faf',
            'gnomad_exomes_faf',
        ],
    )

    # keep previously aggregated transcript consequence field https://github.com/populationgenomics/cpg-flow-seqr-loader/blob/main/src/cpg_seqr_loader/scripts/annotate_cohort.py#L267
    fields_to_keep.append('transcriptConsequenceTerms')

    # keep the max AlphaMissense score for all transcript consequences.
    mt = mt.annotate_rows(
        am_max_score=hl.or_missing(
            hl.is_defined(mt.vep.transcript_consequences),
            hl.max(
                mt.vep.transcript_consequences.filter(lambda tc: hl.is_defined(tc)).map(
                    lambda tc: tc.am_pathogenicity,
                ),
            ),
        ),
    )
    fields_to_keep.append('am_max_score')

    # UTR annotations will be tricky as they are applied per transcript, and you may want to keep both the transcript and result?
    # if you want to keep just the list of 5'UTR predicted consequences you can use some aggregation logic similar to above
    mt = mt.annotate_rows(
        utr_5_consequences=hl.or_missing(
            hl.is_defined(mt.vep.transcript_consequences),
            hl.set(
                mt.vep.transcript_consequences.filter(lambda tc: hl.is_defined(tc))
                .filter(lambda tc: hl.is_defined(tc['5utr_consequence']))
                .map(lambda tc: tc['5utr_consequence']),
            ),
        ),
    )
    fields_to_keep.append('utr_5_consequences')

    # keep AVI scores
    fields_to_keep.append('avis')

    # remove all entries (genotypes) from the dataset where the sample was WT/HomRef
    mt = mt.filter_entries(mt.GT.is_hom_ref(), keep=False)

    # remove variants with low GQ
    min_gq = 20
    mt = mt.filter_entries(min_gq < mt.GQ)

    # remove variants that didn't pass VQSR filters. Variants that 'PASS' are left blank.
    mt = mt.filter_rows(hl.len(mt.filters) == 0)

    # remove common variants
    max_af = 0.01
    mt = mt.filter_rows(max_af > mt.gnomad_genomes.FAF_AF)
    mt = mt.filter_rows(max_af > mt.gnomad_exomes.FAF_AF)

    # remove variants with high AC
    max_ac = 50
    mt = mt.filter_rows(hl.or_missing(hl.len(mt.info.AC) > 0, mt.info.AC[0]) < max_ac)
    # aggregate all sample IDs remaining
    # this would create a new field, `var_samples`, which is a set of all CPG IDs with variants fitting above criteria
    mt = mt.annotate_rows(
        var_samples=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect_as_set(mt.s)),
    )
    fields_to_keep.append('var_samples')

    # adding checkpoint for stability
    checkpoint_path = output_path.replace('.ht', '_checkpoint.mt')
    mt = mt.checkpoint(checkpoint_path, overwrite=True)

    # once we've pulled out all the entry data we want, we can drop all the entries, and all the fields we no longer need
    ht = mt.rows()

    ht = ht.select(*fields_to_keep)

    # write this skinny HT to a new location
    ht.show(5)  # test to make sure logic works
    ht.write(output_path, overwrite=True)

    # print the number of rows
    written_ht = hl.read_table(output_path)
    print(f"Final table size: {written_ht.count()} rows")

    # export as a tsv too?
    tsv_path = output_path.replace('.ht', '.tsv')
    written_ht.export(tsv_path, delimiter='\t')


if __name__ == "__main__":

    # load command line parameters - in and out
    parser = ArgumentParser()
    parser.add_argument("--mt", help="Matrix table to reformat")
    parser.add_argument("--output", help="Write path, format TBD")
    args = parser.parse_args()

    main(args.mt, args.output)
