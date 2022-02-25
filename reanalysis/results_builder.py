"""
All the content relating to presenting the results
"""

from typing import Dict, Set

import pandas as pd

from reanalysis.utils import ReportedVariant, string_format_variant

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
    ):
        """

        :param results_dict:
        :param seqr_lookup:
        :param panelapp_data:
        """
        self.results = results_dict
        self.seqr = seqr_lookup
        self.panelapp = panelapp_data

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
                        # 'csq': variant.var_data.var.INFO.get('csq'),
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
