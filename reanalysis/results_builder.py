"""
All the content relating to presenting the results
"""

from typing import Dict, Set

import pandas as pd

from reanalysis.utils import PedPerson, ReportedVariant, string_format_variant

SEQR_TEMPLATE = (
    '<a href="https://seqr.populationgenomics.org.au/variant_search/'
    'variant/{variant}/family/{family}" target="_blank">{variant}</a>'
)

GNOMAD_TEMPLATE = (
    '<a href="https://gnomad.broadinstitute.org/variant/'
    '{variant}?dataset=gnomad_r3" target="_blank">{value:.5f}</a>'
)
PANELAPP_TEMPLATE = (
    '<a href="https://panelapp.agha.umccr.org/panels/137/gene/{symbol}/"'
    ' target="_blank">{symbol}</a>'
)

ClassColours = {
    1: '<span style="color: #C70039">Class1</span>',
    2: '<span style="color: #FF5733">Class2</span>',
    3: '<span style="color: #FFA533">Class3</span>',
    4: '<span style="color: #33E8FF">Class4</span>',
}


class HTMLBuilder:
    """
    takes the input, makes the output
    """

    def __init__(
        self,
        results_dict: Dict[str, Dict[str, ReportedVariant]],
        seqr_lookup: Dict[str, str],
        panelapp_data: Dict[str, Dict[str, str]],
        csq_string: str,
        pedigree: Dict[str, PedPerson],
    ):  # pylint: disable=too-many-arguments
        """

        :param results_dict:
        :param seqr_lookup:
        :param panelapp_data:
        :param csq_string:
        :param pedigree:
        """
        self.results = results_dict
        self.seqr = seqr_lookup
        self.panelapp = panelapp_data
        self.csq_entries = csq_string.split('|')
        self.pedigree = pedigree

    def get_summary_stats(self) -> str:  # pylint: disable=too-many-locals
        """
        run the numbers across all variant categories
        """

        class_1 = []
        class_2 = []
        class_3 = []
        all_var_strings = set()
        class_1_strings = set()
        class_2_strings = set()
        class_3_strings = set()
        total_count = []
        num_samples = 0
        num_with_variants = 0

        samples_with_no_variants = []

        for sample, entity in self.pedigree.items():
            if not entity.affected:
                continue
            num_samples += 1
            if sample not in self.results.keys():
                samples_with_no_variants.append(sample)
                class_1.append(0)
                class_2.append(0)
                class_3.append(0)
                total_count.append(0)
                continue

            num_with_variants += 1

            # get variants for this sample
            sample_variants = self.results.get(sample)

            # how many variants were attached to this sample?
            total_count.append(len(sample_variants))

            sample_c_1 = 0
            sample_c_2 = 0
            sample_c_3 = 0

            # iterate over the variants
            for variant in sample_variants.values():

                # find all classes associated with this variant
                variant_ints = variant.var_data.get_class_ints()

                # get the string representation
                var_string = string_format_variant(variant.var_data.var)

                # update the set of all unique variants
                all_var_strings.add(var_string)

                # for each class, add to corresponding list and set
                if 1 in variant_ints:
                    sample_c_1 += 1
                    class_1_strings.add(var_string)
                if 2 in variant_ints:
                    sample_c_2 += 1
                    class_2_strings.add(var_string)
                if 3 in variant_ints:
                    sample_c_3 += 1
                    class_3_strings.add(var_string)

            # update the global lists with per-sample counts
            class_1.append(sample_c_1)
            class_2.append(sample_c_2)
            class_3.append(sample_c_3)

        summary_dicts = []

        # turn results into a list of dicts
        for title, deets, unique in [
            ('Total', total_count, all_var_strings),
            ('Class1', class_1, class_1_strings),
            ('Class2', class_2, class_2_strings),
            ('Class3', class_3, class_3_strings),
        ]:
            summary_dicts.append(
                {
                    'Category': title,
                    'Total': sum(deets),
                    'Unique': len(unique),
                    'Peak #/sample': max(deets),
                    'Mean/sample': sum(deets) / len(deets),
                }
            )

        return (
            pd.DataFrame(summary_dicts).to_html(index=False, escape=False),
            samples_with_no_variants,
        )

    def get_csq_from_variant(self, variant: ReportedVariant) -> str:
        """
        populates a single string of all relevant consequences
        """

        csq_set = set()
        csq_info = variant.var_data.var.INFO.get('CSQ')

        for each_csq in csq_info.split(','):
            csq_dict = dict(zip(self.csq_entries, each_csq.split('|')))
            csq_set.update(set(csq_dict['Consequence'].split('&')))
        return ', '.join(csq_set)

    def create_html_tables(self):
        """
        make a table describing the outputs
        :return:
        """
        candidate_dictionaries = {}
        sample_tables = {}

        class_2_genes = set()

        for sample, variants in self.results.items():
            for variant in variants.values():
                variant_class_ints = variant.var_data.get_class_ints()

                if 2 in variant_class_ints:
                    class_2_genes.add(variant.gene)

                var_string = string_format_variant(variant.var_data.var)
                candidate_dictionaries.setdefault(variant.sample, []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'classes': ', '.join(
                            list(
                                map(
                                    lambda x: ClassColours[x],
                                    variant_class_ints,
                                )
                            )
                        ),
                        'symbol': self.panelapp.get(variant.gene).get('symbol'),
                        # 'gene': variant.gene,
                        'csq': self.get_csq_from_variant(variant),
                        'gnomad': GNOMAD_TEMPLATE.format(
                            variant=var_string,
                            value=float(variant.var_data.var.INFO.get('gnomad_af')),
                        ),
                        'MOIs': ','.join(variant.reasons),
                        'support': self.make_seqr_link(
                            var_string=string_format_variant(variant.support_var.var),
                            sample=sample,
                        )
                        if variant.supported
                        else 'N/A',
                    }
                )

        for sample, variants in candidate_dictionaries.items():
            sample_tables[sample] = pd.DataFrame(variants).to_html(
                index=False, render_links=True, escape=False
            )

        return sample_tables, class_2_genes

    def class_2_table(self, class_2_variants: Set[str]) -> str:
        """
        takes all class 2 variants, and documents the panel changes
        this is either a new entity, or an altered MOI
        :param class_2_variants:
        :return:
        """
        current_key = (
            f'MOI in v{self.panelapp["panel_metadata"].get("current_version")}'
        )
        previous_key = (
            f'MOI in v{self.panelapp["panel_metadata"].get("previous_version")}'
        )

        # if we don't have version differences, don't do anything
        if previous_key is None:
            return ''

        gene_dicts = []
        for gene in class_2_variants:
            gene_data = self.panelapp.get(gene)
            gene_dicts.append(
                {
                    'gene': gene,
                    'symbol': PANELAPP_TEMPLATE.format(symbol=gene_data.get('symbol')),
                    previous_key: 'Gene Not Present'
                    if gene_data.get('new')
                    else gene_data.get('old_moi'),
                    current_key: gene_data.get('moi'),
                }
            )
        return pd.DataFrame(gene_dicts).to_html(
            index=False, render_links=True, escape=False
        )

    def make_seqr_link(self, var_string: str, sample: str) -> str:
        """
        either return just the variant as a string, or a seqr link if possible
        :param sample:
        :param var_string:
        :return:
        """
        if sample not in self.seqr:
            return var_string
        return SEQR_TEMPLATE.format(
            variant=var_string,
            family=self.seqr.get(sample),
        )
