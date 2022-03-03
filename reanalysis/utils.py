"""
a collection of classes and methods
which may be shared across reanalysis components
"""


from typing import Any, Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
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


def get_simple_moi(moi: str) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a reduced
    range of cases which can be easily implemented in RD analysis

    Note: Maybe don't do this. Leave the raw data in until the final analysis stage
    :param moi: full PanelApp string
    :return: a simplified representation
    """

    # default to considering both
    # NOTE! The Y chromosome genes all have 'Unknown'!
    panel_app_moi = 'Mono_And_Biallelic'
    if moi is None:
        # exit iteration, return both (all considered)
        return panel_app_moi

    # ideal for match-case, coming to a python 3.10 near you!
    if moi.startswith('BIALLELIC'):
        panel_app_moi = 'Biallelic'
    if moi.startswith("BOTH"):
        panel_app_moi = 'Mono_And_Biallelic'
    if moi.startswith('MONO'):
        panel_app_moi = 'Monoallelic'
    if moi.startswith('X-LINKED'):
        if 'biallelic' in moi:
            panel_app_moi = 'Hemi_Bi_In_Female'
        else:
            panel_app_moi = 'Hemi_Mono_In_Female'

    return panel_app_moi


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
    a prototype already exists in the prototype folder
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


def extract_info(variant: Variant, config: Dict[str, Any]):
    """
    creates an INFO dict by pulling content from the variant info
    keeps a list of dictionaries for each transcript_consequence
    :param variant:
    :param config:
    :return:
    """

    # grab the basic information from INFO
    info_dict = {x: y for x, y in variant.INFO if x in config.get("var_info_keep", [])}

    csq_contents = variant.INFO.get(config.get('csq_field'))
    csq_string = config.get('csq_string')

    # break the mono-CSQ-string into the component headings
    csq_categories = list(map(str.lower, csq_string.split('|')))

    # iterate over all consequences, and make each into a dict
    info_dict['transcript_consequences'] = [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]
    return info_dict


class AnalysisVariant:
    """
    create a variant superclass
    pull all content out of the cyvcf2 object
    we could have a separate implementation for pyvcf,
    or for a direct parser...
    whatever
    """

    def __init__(self, var: Variant, samples: List[str], config: Dict[str, Any]):

        # store the string representation
        self.string = string_format_variant(var)

        # variant coordinate details
        self.chrom = var.CHROM
        self.pos = var.POS
        self.ref = var.REF
        self.alt = var.ALT[0]

        # set the class attributes
        self.class_1 = var.INFO.get('Class1') == 1
        self.class_2 = var.INFO.get('Class2') == 1
        self.class_3 = var.INFO.get('Class3') == 1
        self.class_4 = var.INFO.get('Class4') == 1

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

        self.info = extract_info(variant=var, config=config)

    @property
    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        :return:
        """
        return any([self.class_1, self.class_2, self.class_3, self.class_4])

    @property
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

    @property
    def class_ints(self) -> List[int]:
        """
        get a list of ints representing the classes present on this variant
        for each numerical class, append that number if the class is present
        """
        return [
            integer
            for integer, class_bool in enumerate(
                [
                    self.class_1,
                    self.class_2,
                    self.class_3,
                    self.class_4,
                ],
                1,
            )
            if class_bool
        ]


@dataclass
class ReportedVariant:
    """
    minimal model representing variant categorisation event
    the initial variant
    the MOI passed
    the support (if any)
    """

    sample: str
    gene: str
    var_data: AnalysisVariant
    reasons: Set[str]
    supported: bool
    support_var: Optional[AnalysisVariant] = None
