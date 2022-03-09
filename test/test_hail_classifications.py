"""
unit testing collection for the hail MT methods

aim - reconstruct some fields into a struct
simulate the
"""

import os
import hail as hl

from hail.utils.java import FatalError
import pytest
import pandas as pd

from reanalysis.hail_filter_and_classify import annotate_class_1


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')

"""
try a parametrize with writing then reading an annotation tsv 
"""
titles = [
    'chr',
    'pos',
    'AC',
    'AN',
    'consequence_terms',
    'geneIds',
    'clinvar_sig',
    'clinvar_stars',
    'cadd',
    'revel',
    'lof',
    'sift_score',
    'polyphen_score',
    'mutationtaster',
    'gerp_rs',
    'eigen_phred',
    'gnomad_af',
    'exac_af',
    'expected',
]


class1_keys = ['locus', 'clinvar_sig', 'clinvar_stars']


@pytest.fixture(scope='session')
def hail_matrix():
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - already initiated as GRCh37?')
    return hl.import_vcf(HAIL_VCF, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'values,classified',
    [
        (
            [
                hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'),
                "pathogenic",
                0,
            ],
            False,
        ),
        (
            [
                hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'),
                "pathogenic",
                1,
            ],
            True,
        ),
        (
            [
                hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'),
                "conflicting_interpretations_of_pathogenicity",
                1,
            ],
            False,
        ),
        (
            [
                hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'),
                "benign",
                1,
            ],
            False,
        ),
        (
            [
                hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'),
                "pathogenic&something&else",
                2,
            ],
            True,
        ),
    ],
)
def test_class_1_assignment(values, classified, hail_matrix):
    # cast the input as a dictionary
    row_dict = dict(zip(class1_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            clinvar_stars=hail_table[hail_matrix.locus].clinvar_stars,
        )
    )

    anno_matrix = annotate_class_1(anno_matrix)
    assert anno_matrix.info.Class1.collect() == [classified]
