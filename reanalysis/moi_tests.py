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
from typing import Any, Dict, List, Optional, Tuple

from reanalysis.utils import AnalysisVariant, PedPerson, ReportedVariant


# config keys to use for dominant MOI tests
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom', 'exac_ac_hom'}


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

    # search for a partner-hit using the string-format of this variant
    if first_variant.string in comp_hets.get(sample):
        paired_variant = comp_hets[sample][first_variant.string]
        # check that the class combination is supported
        if all([first_variant.class_4_only, paired_variant.class_4_only]):
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
        self,
        pedigree: Dict[str, PedPerson],
        target_moi: str,
        config: Dict[str, Any],
    ):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant

        :param pedigree:
        :param target_moi:
        :param config:
        """

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi in ['Monoallelic']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree, config=config),
            ]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree, config=config),
                RecessiveAutosomal(pedigree=pedigree, config=config),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [RecessiveAutosomal(pedigree=pedigree, config=config)]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [
                XRecessive(pedigree=pedigree, config=config),
                XDominant(pedigree=pedigree, config=config),
            ]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [XRecessive(pedigree=pedigree, config=config)]

        elif target_moi == 'Y_Chrom_Variant':
            self.filter_list = [YHemi(pedigree=pedigree, config=config)]

        else:
            raise Exception(f'MOI type {target_moi} is not addressed in MOI')

    def run(
        self, principal_var, comp_hets: Dict[str, Dict[str, AnalysisVariant]], ensg: str
    ) -> List[ReportedVariant]:
        """
        run method - triggers each relevant inheritance model
        :param principal_var: the variant we are focused on
        :param comp_hets: all compound hets
        :param ensg: the specific gene ID we're looking at
            this can be pulled from the annotations, but could be required
            if we stop splitting out each consequence into a new VCF row
        :return:
        """
        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(model.run(principal_var, comp_hets, ensg))
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

    def __init__(
        self, pedigree: Dict[str, PedPerson], config: Dict[str, Any], applied_moi: str
    ):
        """
        no values to establish?
        maybe constants in the base class
        """
        if applied_moi is None:
            raise Exception("An applied MOI needs to reach the Base Class")
        self.pedigree = pedigree
        self.config = config
        self.applied_moi = applied_moi

    @abstractmethod
    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        :param principal_var:
        :param comp_hets:
        :param ensg:
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
        config: Dict[str, Any],
        applied_moi: str = 'Autosomal Dominant',
    ):
        """

        :param pedigree: not yet implemented
        :param config:
        :param applied_moi:
        """
        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        self.hom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        super().__init__(pedigree=pedigree, config=config, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """
        simplest
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param comp_hets:
        :param ensg:
        :return:
        """

        classifications = []

        # more stringent Pop.Freq checks for dominant
        if (
            principal_var.info.get('gnomad_af') >= self.ad_threshold
            or any(
                [
                    principal_var.info.get(hom_key) >= self.hom_threshold
                    for hom_key in INFO_HOMS
                ]
            )
            or principal_var.info.get('gnomad_ac') >= self.ac_threshold
        ):
            return classifications

        # autosomal dominant doesn't require support
        for sample_id in [
            sam for sam in principal_var.het_samples if self.is_affected(sam)
        ]:
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=ensg,
                    var_data=principal_var,
                    reasons={self.applied_moi},
                    supported=False,
                )
            )

        return classifications


class RecessiveAutosomal(BaseMoi):
    """
    pass
    """

    def __init__(
        self,
        pedigree: Dict[str, PedPerson],
        config: Dict[str, Any],
        applied_moi: str = 'Autosomal Recessive',
    ):
        """ """
        self.hom_threshold = config.get(GNOMAD_REC_HOM_THRESHOLD)
        super().__init__(pedigree=pedigree, config=config, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """
        valid if present as hom, or compound het
        counts as being phased if a compound het is split between parents
        Clarify if we want to consider a homozygous variant as 2 hets
        :param principal_var:
        :param comp_hets:
        :param ensg:
        :return:
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if any(
            [
                principal_var.info.get(hom_key) >= self.hom_threshold
                for hom_key in INFO_HOMS
            ]
        ):
            return classifications

        # homozygous is relevant directly
        for sample_id in [
            sam for sam in principal_var.hom_samples if self.is_affected(sam)
        ]:

            # # check for either parent being unaffected - skip this sample
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=ensg,
                    var_data=principal_var,
                    reasons={f'{self.applied_moi} Homozygous'},
                    supported=False,
                )
            )

        # if hets are present, try and find support
        for sample_id in [
            sam for sam in principal_var.het_samples if self.is_affected(sam)
        ]:

            passes, partner = check_for_second_hit(principal_var, comp_hets, sample_id)
            if passes:
                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        gene=ensg,
                        var_data=principal_var,
                        reasons={f'{self.applied_moi} Compound-Het'},
                        supported=True,
                        support_var=partner,
                    )
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

    def __init__(
        self,
        pedigree: Dict[str, PedPerson],
        config: Dict[str, Any],
        applied_moi: str = 'X_Dominant',
    ):
        """
        accept male hets and homs, and female hets without support
        :param pedigree:
        :param config:
        :param applied_moi:
        """
        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        self.hom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        super().__init__(pedigree=pedigree, config=config, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """
        if variant is present and sufficiently rare, we take it

        :param principal_var:
        :param comp_hets:
        :param ensg:
        :return:
        """
        classifications = []

        if principal_var.chrom.replace('chr', '').lower() != 'x':
            logging.error(
                'X-Chromosome MOI given for variant on %s', principal_var.chrom
            )

        # more stringent Pop.Freq checks for dominant
        if (
            principal_var.info.get('gnomad_af') >= self.ad_threshold
            or any(
                [
                    principal_var.info.get(hom_key) >= self.hom_threshold
                    for hom_key in INFO_HOMS
                ]
            )
            or principal_var.info.get('gnomad_ac') >= self.ac_threshold
        ):
            return classifications

        # due to legacy caller issues, we expect wrongly called HOM Males
        # X-relevant, we separate out male and females
        males = [
            sam
            for sam in principal_var.het_samples.union(principal_var.hom_samples)
            if self.pedigree.get(sam).male
        ]

        het_females = [
            sam for sam in principal_var.het_samples if not self.pedigree.get(sam).male
        ]

        # dominant doesn't care about support
        # for all the samples which are affected
        # (assumption that probands are affected)
        for sample_id in males + het_females:
            if self.pedigree.get(sample_id).male:
                reason = f'{self.applied_moi} Male'
            else:
                reason = f'{self.applied_moi} Female'
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=ensg,
                    var_data=principal_var,
                    reasons={reason},
                    supported=False,
                )
            )
        return classifications


class XRecessive(BaseMoi):
    """
    for males accept het** - male variants HOM because GATK
    effectively the same as AutosomalDominant?
    """

    def __init__(
        self,
        pedigree: Dict[str, PedPerson],
        config: Dict[str, Any],
        applied_moi: str = 'X_Recessive',
    ):
        """ """

        self.hom_dom_threshold = config.get(GNOMAD_DOM_HOM_THRESHOLD)
        self.hom_rec_threshold = config.get(GNOMAD_REC_HOM_THRESHOLD)

        super().__init__(pedigree=pedigree, config=config, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """

        :param principal_var:
        :param comp_hets:
        :param ensg:
        :return:
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if any(
            [
                principal_var.info.get(hom_key) >= self.hom_dom_threshold
                for hom_key in INFO_HOMS
            ]
        ):
            return classifications

        # X-relevant, we separate out male and females
        males = [
            sam
            for sam in principal_var.het_samples.union(principal_var.hom_samples)
            if self.pedigree.get(sam).male
        ]

        # get female calls in 2 categories
        het_females = [
            sam for sam in principal_var.het_samples if not self.pedigree.get(sam).male
        ]

        hom_females = [
            sam for sam in principal_var.hom_samples if not self.pedigree.get(sam).male
        ]

        # if het females are present, try and find support
        for sample_id in het_females:

            passes, partner = check_for_second_hit(principal_var, comp_hets, sample_id)
            if passes:
                classifications.append(
                    ReportedVariant(
                        sample=sample_id,
                        gene=ensg,
                        var_data=principal_var,
                        reasons={f'{self.applied_moi} Compound-Het Female'},
                        supported=True,
                        support_var=partner,
                    )
                )

        # remove from analysis if too many homs are present in population databases
        if any(
            [
                principal_var.info.get(hom_key) >= self.hom_rec_threshold
                for hom_key in INFO_HOMS
            ]
        ):
            return classifications

        # find all het males and hom females
        # assumption that the sample can only be hom if female?
        for sample_id in males + hom_females:
            if self.pedigree.get(sample_id).male:
                reason = f'{self.applied_moi} Male'
            else:
                reason = f'{self.applied_moi} Female'
            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=ensg,
                    var_data=principal_var,
                    reasons={reason},
                    supported=False,
                )
            )
        return classifications


class YHemi(BaseMoi):
    """
    so... we should flag any female calls?
    otherwise treat as AD
    not expecting to use this
    """

    def __init__(
        self,
        pedigree: Dict[str, PedPerson],
        config: Dict[str, Any],
        applied_moi: str = 'Y_Hemi',
    ):
        """

        :param pedigree:
        :param config:
        :param applied_moi:
        """

        self.ad_threshold = config.get(GNOMAD_RARE_THRESHOLD)
        self.ac_threshold = config.get(GNOMAD_AD_AC_THRESHOLD)
        super().__init__(pedigree=pedigree, config=config, applied_moi=applied_moi)

    def run(
        self,
        principal_var: AnalysisVariant,
        comp_hets: Dict[str, Dict[str, AnalysisVariant]],
        ensg: str,
    ) -> List[ReportedVariant]:
        """
        flag calls on Y which are Hom (maybe ok?) or female (bit weird)
        :param principal_var:
        :param comp_hets:
        :param ensg:
        :return:
        """
        classifications = []

        # more stringent Pop.Freq checks for dominant
        if (
            principal_var.info.get('gnomad_af') >= self.ad_threshold
            or principal_var.info.get('gnomad_ac') >= self.ac_threshold
        ):
            return classifications

        # y chrom... called as hom? Shouldn't be possible
        # half expecting this with GATK...
        for sample_id in principal_var.hom_samples:
            logging.warning('Sample %s is a hom call on Y', sample_id)

        # we don't expect any confident Y calls in females
        for sample_id in principal_var.het_samples.union(principal_var.hom_samples):
            if not self.pedigree.get(sample_id).male:
                logging.error('Sample %s is a female with call on Y', sample_id)

            classifications.append(
                ReportedVariant(
                    sample=sample_id,
                    gene=ensg,
                    var_data=principal_var,
                    reasons={self.applied_moi},
                    supported=False,
                )
            )

        return classifications
