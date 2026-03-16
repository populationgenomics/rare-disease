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
                mt.vep.transcript_consequences
                .filter(lambda tc: hl.is_defined(tc))
                .map(lambda tc: tc.am_pathogenicity)
            )
        )
    )
    fields_to_keep.append('am_max_score')
    
    # UTR annotations will be tricky as they are applied per transcript, and you may want to keep both the transcript and result?
    # if you want to keep just the list of 5'UTR predicted consequences you can use some aggregation logic similar to above
    mt = mt.annotate_rows(
        utr_5_consequences=hl.or_missing(
            hl.is_defined(mt.vep.transcript_consequences),
            hl.set(
                mt.vep.transcript_consequences
                .filter(lambda tc: hl.is_defined(tc))
                .filter(lambda tc: hl.is_defined(tc['5utr_consequence']))
                .map(lambda tc: tc['5utr_consequence'])
            )
        )
    )
    fields_to_keep.append('utr_5_consequences')

	# keep AVI scores
    fields_to_keep.append('avis')

    # remove all entries (genotypes) from the dataset where the sample was WT/HomRef
    mt = mt.filter_entries(mt.GT.is_hom_ref(), keep=False)

    # remove variants with low GQ
    mt = mt.filter_entries(mt.GQ > 20)

    # remove variants that didn't pass VQSR filters. Variants that 'PASS' are left blank.
    mt = mt.filter_rows(hl.len(mt.filters) == 0)

    # remove common variants
    mt = mt.filter_rows(mt.gnomad_genomes.FAF_AF < 0.1)
    mt = mt.filter_rows(mt.gnomad_exomes.FAF_AF < 0.1)

	# remove variants with high AC
    mt = mt.filter_rows(
        hl.or_missing(
            hl.len(mt.info.AC) > 0, 
            mt.info.AC[0]
        ) < 50
    )
    # aggregate all sample IDs remaining
    # this would create a new field, `var_samples`, which is a set of all CPG IDs with variants fitting above criteria
    mt = mt.annotate_rows(
        var_samples=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect_as_set(mt.s))
    )
    fields_to_keep.append('var_samples')

    # it probably makes sense to keep all genotypes fitting your strict criteria here
    # later when you want to make specific choices, such as 'only affected', or 'only with RNA data'
    # you can add a list of samples, make it into a hail object, and do some matches, e.g.
    # rna_samples = ['cpg1', 'cpg3', ]
    ## turn the variable into an expression (hail likes expressions...)
    # hl_samples = hl.literal(rna_samples)

    ## keep all rows where there's a sample with matched RNA
    # ht.filter(hl.len(hl_samples.intersection(ht.var_samples)) > 0)

    ## remove rows unless they're only related to affected participants
    # affected_sams = [...]
    # ht.filter(ht.var_samples.is_subset(affected_sams))

    # there's much more comprehensive logic for entries here: https://github.com/populationgenomics/cpg-flow-seqr-loader/blob/main/src/cpg_seqr_loader/scripts/annotate_dataset.py#L13
    # annotate_dataset generates a number of columns which retains more sample level information, e.g. sample GQ -> samples_gq.5_to_10, samples_gq.10_to_15, samples_gq.15_to_20
    # this might be useful, it might be overkill - you might want to duplicate some of that into this script
    
	# adding checkpoint for stability
    checkpoint_path = output_path.replace('.ht', '_checkpoint.mt')
    mt = mt.checkpoint(checkpoint_path, overwrite=True)

    # once we've pulled out all the entry data we want, we can drop all the entries, and all the fields we no longer need
    ht = mt.rows()

    ht = ht.select(*fields_to_keep)

    # write this skinny HT to a new location
    ht.show(5) #test to make sure logic works
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
