"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Conceptually borrowing from the Genomics England tiering approach

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant,
e.g. if the MOI is Dominant, we may also be interested in Recessive (?)
e.g. if the MOI is X-linked Dom, we also search X-linked Recessive (?)
    - relevant for females, if the pedigree is loaded
    - also depends on the accuracy of a male hemi call

This will come down to clinician preference, some may exclusively
prefer dominant model = monogenic search

This model does not apply anything to MT, I expect those to default
to a Monogenic MOI
"""


import logging
from abc import abstractmethod
from typing import Dict, List, Optional, Tuple


from reanalysis.utils import (
    AnalysisVariant,
    PedPerson,
    c4_only,
    string_format_variant,
)


def check_for_second_hit(
    first_variant: AnalysisVariant,
    comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    sample: str,
) -> Tuple[bool, Optional[AnalysisVariant]]:
    """
    checks for a second hit to this compound het
    :param first_variant:
    :param comp_hets:
    :param sample:
    :return:
    """

    response = False, None
    if sample not in comp_hets.keys():
        return response

    # find the string repr of this variant
    variant_string = string_format_variant(first_variant.var)

    # search for a partner-hit
    if variant_string in comp_hets.get(sample):
        paired_variant = comp_hets[sample][variant_string]
        # check that the class combination is supported
        if all([c4_only(first_variant.var), c4_only(paired_variant.var)]):
            logging.info(
                'Variant pairing discarded as both sides of comp-het were C4 only'
            )
        else:
            response = True, paired_variant
    return response


class MOIRunner:
    """
    pass
    """

    def __init__(
        self, ad_threshold: float, pedigree: Dict[str, PedPerson], target_moi: str
    ):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant

        :param ad_threshold: more stringent AF threshold for dominant inheritance
        :param pedigree:
        :param target_moi:
        """

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi in ['Monoallelic']:
            self.filter_list = [
                DominantAutosomal(pedigree, ad_threshold),
            ]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree, ad_threshold),
                RecessiveAutosomal(pedigree),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [RecessiveAutosomal(pedigree)]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [XRecessive(pedigree), XDominant(pedigree, ad_threshold)]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [XRecessive(pedigree)]

        elif target_moi == 'Y_Chrom_Variant':
            self.filter_list = [YHemi(pedigree, ad_threshold)]

        else:
            raise Exception(f'MOI type {target_moi} is not addressed in MOI')

    def run(self, principal_var, comp_hets: Dict[str, Dict[str, AnalysisVariant]]):
        """
        run method - triggers each relevant inheritance model
        :param principal_var: the variant we are focused on
        :param comp_hets: all compound hets
        :return:
        """
        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(model.run(principal_var, comp_hets))
        return moi_matched

    def send_it(self):
        """
        to stop pylint complaining
        :return:
        """
        print(f'Yaaas {self.filter_list}')


class BaseMoi:
    """
    pass
    """

    def __init__(self, pedigree: Dict[str, PedPerson], applied_moi: str):
        """
        no values to establish?
        maybe constants in the base class
        """
        if applied_moi is None:
            raise Exception("An applied MOI needs to reach the Base Class")
        self.applied_moi = applied_moi
        self.pedigree = pedigree

    @abstractmethod
    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> Tuple[str, str, Tuple[AnalysisVariant]]:
        """
        run over all the applicable inheritance patterns and finds good fits
        :return:
        """

    def is_affected(self, sample_id: str) -> bool:
        """
        take a sample ID and check if they are affected
        :param sample_id:
        :return:
        """
        if self.pedigree.get(sample_id).affected:
            return True
        return False

    # not implemented for singleton analysis
    # def get_affected_parents(self, sample_id: str) -> Set[str]:
    #     """
    #     for this sample, find any affected parents
    #     port this back into the Participant object?!
    #     prevent re-computation... but really only make that more complex if necessary
    #     :param sample_id:
    #     :return:
    #     """
    #     participant = self.pedigree.get(sample_id)
    #
    #     return {
    #         parent.details.sample_id
    #         for parent in [participant.mother, participant.father]
    #         if parent is not None and parent.details.affected
    #     }


class DominantAutosomal(BaseMoi):
    """
    This class can also be called by the X-linked Dominant, in which case the
    Applied_MOI by name is overridden
    """

    def __init__(
        self,
        pedigree: Dict[str, PedPerson],
        ad_threshold: float,
        applied_moi: str = 'Autosomal_Dominant',
    ):
        """

        :param pedigree: not yet implemented
        :param ad_threshold:
        :param applied_moi:
        """
        self.ad_thresh = ad_threshold
        super().__init__(pedigree, applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> List[Tuple[str, str, Tuple[AnalysisVariant]]]:
        """
        simplest
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param comp_hets:
        :return:
        """

        classifications = []

        # apply a more stringent AF threshold for dominant
        if principal_var.var.INFO.get('gnomad_af') >= self.ad_thresh:
            return classifications

        # autosomal dominant doesn't care about support
        # for all the samples which are affected
        # assumption that probands are affected
        for sample_id in [
            sam for sam in principal_var.het_samples if self.is_affected(sam)
        ]:
            classifications.append((sample_id, 'Autosomal Dominant', (principal_var,)))

        return classifications


class RecessiveAutosomal(BaseMoi):
    """
    pass
    """

    def __init__(self, pedigree, applied_moi: str = 'Autosomal_Recessive'):
        """ """
        super().__init__(pedigree, applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> List[Tuple[str, str, Tuple[AnalysisVariant]]]:
        """
        valid if present as hom, or comound het
        counts as being phased if a compound het is split between parents
        Clarify if we want to consider a homozygous variant as 2 hets
        :param principal_var:
        :param comp_hets:
        :return:
        """

        classifications = []

        # homozygous is relevant directly
        for sample_id in [
            sam for sam in principal_var.hom_samples if self.is_affected(sam)
        ]:

            # # check for either parent being unaffected - skip this sample
            classifications.append(
                (sample_id, 'Autosomal Homozygous', (principal_var,))
            )

        # if hets are present, try and find support
        for sample_id in [
            sam for sam in principal_var.het_samples if self.is_affected(sam)
        ]:

            passes, partner = check_for_second_hit(principal_var, comp_hets, sample_id)
            if passes:
                logging.info('comp-het found: %s', {repr(partner)})
                classifications.append(
                    (sample_id, 'Autosomal Compound-Het', (principal_var, partner))
                )

        return classifications


class XDominant(BaseMoi):
    """
    for males and females, accept het
    effectively the same as DominantAutosomal?
    just override the type, and use AD

    GATK MALES ARE CALLED HOM (OR HET pseudo-autosomal?)
    re-implement here, but don't permit Male X-Homs
    """

    def __init__(self, pedigree, ad_threshold: float, applied_moi: str = 'X_Dominant'):
        """
        accept male hets and homs, and female hets without support
        :param pedigree:
        :param applied_moi:
        """
        self.ad_thresh = ad_threshold
        super().__init__(pedigree, applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> List[Tuple[str, str, Tuple[AnalysisVariant]]]:
        """
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param comp_hets:
        :return:
        """
        classifications = []

        if principal_var.var.CHROM.replace('chr', '').lower() != 'x':
            logging.error(
                'X-Chromosome MOI given for variant on %s', principal_var.var.CHROM
            )

        # apply a more stringent AF threshold for dominant
        if principal_var.var.INFO.get('gnomad_af') >= self.ad_thresh:
            return classifications

        # due to legacy caller issues, we expect wrongly called HOM Males
        # X-relevant, we separate out male and females
        males = [
            sam for sam in principal_var.variant_samples if self.pedigree.get(sam).male
        ]

        het_females = [
            sam for sam in principal_var.het_samples if not self.pedigree.get(sam).male
        ]

        # dominant doesn't care about support
        # for all the samples which are affected
        # assumption that probands are affected
        for sample_id in males + het_females:
            if self.pedigree.get(sample_id).male:
                reason = 'X-Hemi Male'
            else:
                reason = 'X-Het Female'
            classifications.append(
                (sample_id, f'{reason}, Autosomal Dominant', (principal_var,))
            )

        return classifications


class XRecessive(BaseMoi):
    """
    for males accept het
    effectively the same as AutosomalDominant?
    """

    def __init__(self, pedigree, applied_moi: str = 'X_Recessive'):
        """ """
        super().__init__(pedigree, applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> List[Tuple[str, str, Tuple[AnalysisVariant]]]:
        """

        :param principal_var:
        :param comp_hets:
        :return:
        """

        classifications = []

        # X-relevant, we separate out male and females
        # GATK calls X-males as Hom
        males = [
            sam for sam in principal_var.variant_samples if self.pedigree.get(sam).male
        ]

        # get female calls in 2 categories
        het_females = [
            sam for sam in principal_var.het_samples if not self.pedigree.get(sam).male
        ]

        hom_females = [
            sam for sam in principal_var.hom_samples if not self.pedigree.get(sam).male
        ]

        # find all het males and hom females
        # assumption that the sample can only be hom if female?
        for sample_id in males + hom_females:
            if self.pedigree.get(sample_id).male:
                reason = 'X-Hemi Male'
            else:
                reason = 'X-Hom Female'
            classifications.append((sample_id, reason, (principal_var,)))

        # if het females are present, try and find support
        for sample_id in het_females:

            passes, partner = check_for_second_hit(principal_var, comp_hets, sample_id)
            if passes:
                classifications.append(
                    (sample_id, 'X-Hom Compound-Het Female', (principal_var, partner))
                )
        return classifications


class YHemi(BaseMoi):
    """
    so... we should flag any female calls?
    otherwise treat as AD
    not expecting to use this
    """

    def __init__(self, pedigree, ad_threshold: float, applied_moi: str = 'Y_Hemi'):
        """

        :param pedigree:
        :param ad_threshold:
        :param applied_moi:
        """
        self.ad_thresh = ad_threshold
        super().__init__(pedigree, applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
    ) -> List[Tuple[str, str, Tuple[AnalysisVariant]]]:
        """
        flag calls on Y which are Hom (maybe ok?) or female (bit weird)
        :param principal_var:
        :param comp_hets:
        :return:
        """
        classifications = []

        # check that the variant is rare enough to be Dominant
        if principal_var.var.INFO.get('gnomad_af') >= self.ad_thresh:
            return classifications

        # y chrom... called as hom? Shouldn't be possible
        # half expecting this with GATK...
        for sample_id in principal_var.hom_samples:
            logging.warning('Sample %s is a hom call on Y', sample_id)

        # we don't expect any confident Y calls in females
        for sample_id in principal_var.variant_samples:
            if not self.pedigree.get(sample_id).male:
                logging.error('Sample %s is a female with call on Y', sample_id)

            classifications.append((sample_id, 'Y-Chrom Het Male', (principal_var,)))

        return classifications
