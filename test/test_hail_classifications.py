"""
unit testing collection for the hail MT methods
"""

import os
import hail as hl


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
CLASS_VCF = os.path.join(INPUT, 'single_var.vcf.bgz')

class_1_loci = [hl.Locus(contig='chr1', position=1045488, reference_genome='GRCh38')]
class_1_loci_b = [hl.Locus(contig='chr1', position=12497580, reference_genome='GRCh38')]
class_1_loci_c = [
    hl.Locus(contig='chr7', position=107236518, reference_genome='GRCh38')
]
class_1_loci_d = [
    hl.Locus(contig='chr19', position=11446352, reference_genome='GRCh38')
]
class_1_loci_e = [hl.Locus(contig='chrX', position=71108336, reference_genome='GRCh38')]
