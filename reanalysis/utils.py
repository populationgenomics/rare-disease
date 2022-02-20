"""
a collection of classes and methods
which may be shared across reanalysis components
"""


from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from enum import Enum
from csv import DictReader

from cyvcf2 import Variant


COMP_HET_VALUES = ['sample', 'gene', 'id', 'chrom', 'pos', 'ref', 'alt']
VARIANT_STRING_TEMPLATE = '{}-{}-{}-{}'
COMP_HET_TEMPLATE = f'{VARIANT_STRING_TEMPLATE}-{{}}'

HOMREF = 0
HETALT = 1
UNKNOWN = 2
HOMALT = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES = {HOMREF, UNKNOWN}


def string_format_variant(var: Variant, transcript: Optional[bool] = False) -> str:
    """
    generates a gnomad and seqr format variant string
    :param var:
    :param transcript: if present, include transcript ID
    :return:
    """

    var_string = VARIANT_STRING_TEMPLATE.format(
        var.CHROM.replace('chr', ''), var.POS, var.REF, var.ALT[0]
    )
    # if transcript was sent, include in the string
    if transcript:
        var_string = f'{var_string}-{var.INFO.get("transcript_id")}'

    return var_string


@dataclass
class PedPerson:
    """
    holds attributes about a single PED file entry
    this will need to be enhanced for family analysis
    """

    sample: str
    male: bool
    affected: bool


def parse_ped_simple(ped: str) -> Dict[str, PedPerson]:
    """
    take individual attributes - sample ID, sex, affected
    :param ped: path to the ped file
    :return:
    """

    ped_dict: Dict[str, PedPerson] = {}
    with open(ped, 'r', encoding='utf-8') as handle:
        for line in DictReader(handle, delimiter="\t"):

            # slot in the sample ID and two Booleans
            ped_dict[line['Individual ID']] = PedPerson(
                line['Individual ID'], line['Sex'] == '1', line['Affected'] == '2'
            )
    return ped_dict


class AnalysisVariant:
    """
    create a variant object with auto-incrementing ID
    """

    def __init__(self, var: Variant, samples: List[str]):
        # full multi-sample variant
        self.var: Variant = var

        # set the class attributes
        self.class_1 = var.INFO.get('Class1') == '1'
        self.class_2 = var.INFO.get('Class2') == '1'
        self.class_3 = var.INFO.get('Class3') == '1'
        self.class_4 = var.INFO.get('Class4') == '1'

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        :return:
        """
        return any([self.class_1, self.class_2, self.class_3, self.class_4])

    def class_4_only(self) -> bool:
        """
        checks that the variant was only class 4
        :return:
        """
        return self.class_4 and not any(
            [
                self.class_1,
                self.class_2,
                self.class_3,
            ]
        )


def get_non_ref_samples(
    variant: Variant, samples: List[str]
) -> Tuple[Set[str], Set[str]]:
    """
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    return het, hom, and the union of het and hom

    maybe something different would be more versatile
    e.g.
    hets = {
        sample: '0/1',
    }
    where the genotype and phased (|) vs unphased (/)
    is determined from the variant.genotypes attribute
    This would make it trivial for the final output to
    have an accurate representation of the parsed GT
    without having to regenerate the string rep.
    :param variant:
    :param samples:
    :return:
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types):

        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == 1:
            het_samples.add(sam)
        if genotype_int == 3:
            hom_samples.add(sam)

    return het_samples, hom_samples


@dataclass
class ReportedVariant:
    """
    an attempt to describe a model variant
    all the details required for a report
    :return:
    """


class SimpleMOI(Enum):
    """
    enumeration to simplify the panelapp MOIs
    """

    BOTH = 'Mono_And_Biallelic'
    BIALLELIC = 'Biallelic'
    MONOALLELIC = 'Monoallelic'
    X_MONO = 'Hemi_Mono_In_Female'
    X_BI = 'Hemi_Bi_In_Female'
    Y_MONO = 'Y_Chrom_Variant'
    UNKNOWN = 'Unknown'


class AppliedMoi(Enum):
    """
    the different inheritance patterns to apply
    assumed complete penetrance only for now
    """

    AUTO_DOM = 'Autosomal_Dominant'
    AUTO_REC = 'Autosomal_Recessive'
    X_DOM = 'X_Dominant'
    X_REC = 'X_Recessive'
    Y_HEMI = 'Y_Hemi'
