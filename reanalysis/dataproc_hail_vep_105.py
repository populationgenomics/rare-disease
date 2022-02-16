#!/usr/bin/env python3


"""
Run VEP on the hail mt
"""


import click
import hail as hl


@click.command()
@click.option('--matrix', required=True, help='Hail matrix table to run VEP on')
@click.option('--output', required=True, help='Hail matrix table to write out')
def main(matrix: str, output: str):
    """
    Run vep using main.py wrapper
    :param matrix: input path
    :param output: output path
    """

    hl.init(default_reference='GRCh38')
    matrix_data = hl.read_matrix_table(matrix)

    # a couple of annotation-independent filtering steps to
    # avoid annotating variants we won't consider anyway

    # hard filter for quality and abundance in the joint call
    # locus restriction to keep this tiny
    matrix_data = matrix_data.filter_rows(
        (matrix_data.info.AC <= 2)
        & (matrix_data.filters.length() == 0)
        & (hl.len(matrix_data.alleles) == 2)
        & (matrix_data.alleles[1] != '*')
        & (matrix_data.locus.contig == 'chr22')
    )

    # throw in a repartition here (annotate even chunks in parallel)
    matrix_data = matrix_data.repartition(20, shuffle=False)

    # now run sexy VEP 105 annotation
    vep = hl.vep(matrix_data, config='file:///vep_data/vep-gcloud.json')
    hl.export_vcf(
        vep,
        output,
        tabix=True,
    )
    # vep.write(output)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
