"""
Methods for taking the final output and generating static report content
"""

from typing import Dict, List, Set, Tuple, Union

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

STRONG_STRING = '<strong>{content}</strong>'
COLOUR_STRING = '<span style="color: {colour}">{content}</span>'


class HTMLBuilder:
    """
    takes the input, makes the output
    """

    def __init__(
        self,
        results_dict: Dict[str, Dict[str, ReportedVariant]],
        seqr_lookup: Dict[str, str],
        panelapp_data: Dict[str, Dict[str, str]],
        config: Dict[str, Union[str, Dict[str, str]]],
        pedigree: Dict[str, PedPerson],
    ):  # pylint: disable=too-many-arguments
        """

        :param results_dict:
        :param seqr_lookup:
        :param panelapp_data:
        :param config:
        :param pedigree:
        """
        self.results = results_dict
        self.seqr = seqr_lookup
        self.panelapp = panelapp_data
        self.config = config
        self.csq_entries = config.get('csq_string').split('|')
        self.pedigree = pedigree
        self.colours = self.set_up_colours()

    def set_up_colours(self):
        """
        generates the dictionary of colour strings
        """

        return {
            key: COLOUR_STRING.format(
                colour=value, content=STRONG_STRING.format(content=f'Class{key}')
            )
            for key, value in self.config["colours"].items()
        }

    def get_summary_stats(
        self,
    ) -> Tuple[str, List[str]]:  # pylint: disable=too-many-locals
        """
        run the numbers across all variant categories
        """

        class_count = {'1': [], '2': [], '3': [], 'all': []}
        class_strings = {'1': set(), '2': set(), '3': set(), 'all': set()}

        samples_with_no_variants = []

        for sample, entity in self.pedigree.items():
            if not entity.affected:
                continue
            if sample not in self.results.keys():
                samples_with_no_variants.append(sample)

                # update all indices; 0 variants for this sample
                for class_list in class_count.values():
                    class_list.append(0)
                continue

            # get variants for this sample
            sample_variants = self.results.get(sample)

            # how many variants were attached to this sample?
            class_count['all'].append(len(sample_variants))

            # create a per-sample object to track variants for each class
            sample_count = {'1': 0, '2': 0, '3': 0}

            # iterate over the variants
            for variant in sample_variants.values():

                # get the string representation
                var_string = string_format_variant(variant.var_data.var)

                # find all classes associated with this variant
                # for each class, add to corresponding list and set
                for class_value in list(map(str, variant.var_data.get_class_ints())):
                    if class_value == '4':
                        continue
                    sample_count[class_value] += 1
                    class_strings[class_value].add(var_string)

                # update the set of all unique variants
                class_strings['all'].add(var_string)

            # update the global lists with per-sample counts
            for key, value in sample_count.items():
                class_count[key].append(value)

        summary_dicts = [
            {
                'Category': title,
                'Total': sum(class_count[key]),
                'Unique': len(class_strings[key]),
                'Peak #/sample': max(class_count[key]),
                'Mean/sample': sum(class_count[key]) / len(class_count[key]),
            }
            for title, key in [
                ('Total', 'all'),
                ('Class1', '1'),
                ('Class2', '2'),
                ('Class3', '3'),
            ]
        ]

        return (
            pd.DataFrame(summary_dicts).to_html(index=False, escape=False),
            samples_with_no_variants,
        )

    def write_html(self, output_path: str):
        """
        uses the class objects to create the HTML tables
        writes all content to the output path
        """

        summary_table, zero_classified_samples = self.get_summary_stats()
        html_tables, class_2_genes = self.create_html_tables()
        class_2_table = self.class_2_table(class_2_genes)

        with open(output_path, 'w', encoding='utf-8') as handle:
            handle.write('<head>\n</head>\n<body>\n')
            handle.write('<h3>MOI changes used for Class 2</h3>')
            handle.write(class_2_table)
            handle.write('<br/>')

            handle.write(
                f'<h3>Samples without Classified Variants ({len(zero_classified_samples)})</h3>'
            )
            if len(zero_classified_samples) > 0:
                handle.write(f'<h5>{", ".join(zero_classified_samples)}</h3>')
            handle.write('<br/>')

            handle.write('<h3>Per-Class summary</h3>')
            handle.write(summary_table)
            handle.write('<br/>')

            handle.write('<h1>Per Sample Results</h1>')
            for sample, table in html_tables.items():
                handle.write(fr'<h3>Sample: {sample}</h3>')
                handle.write(table)
            handle.write('\n</body>')

    def colour_csq(self, all_csq: Set[str], mane_csq: Set[str]) -> str:
        """
        takes the collection of all consequences, and MANE csqs
        if a CSQ occurs on MANE, write in bold,
        if non-MANE, write in red
        return the concatenated string
        """
        csq_strings = []
        for csq in all_csq:
            # bold, in Black
            if csq in mane_csq:
                csq_strings.append(STRONG_STRING.format(content=csq))
            # bold, and red
            else:
                csq_strings.append(
                    self.colours["red_consequence"].replace('Classred_consequence', csq)
                )

        return ", ".join(csq_strings)

    def get_csq_details(self, variant: ReportedVariant) -> Tuple[str, str]:
        """
        populates a single string of all relevant consequences
        UPDATE - take MANE into account
        """

        csq_set = set()
        mane_transcript = set()
        mane_csq = set()
        csq_info = variant.var_data.var.INFO.get('CSQ')

        # iterate over all consequences, special care for MANE
        for each_csq in csq_info.split(','):
            csq_dict = dict(zip(self.csq_entries, each_csq.split('|')))
            row_csq = set(csq_dict['Consequence'].split('&'))

            # record the transcript ID(s), and CSQ(s)
            # we only expect 1, but set operations are useful
            if csq_dict.get('MANE_select') != '':
                mane_csq.update(row_csq)
                mane_transcript.add(csq_dict.get('MANE_select'))
            csq_set.update(row_csq)

        # we only ever expect one... but this would make the addition of plus_clinical easy
        mane_string = STRONG_STRING.format(content=', '.join(mane_csq))

        return self.colour_csq(csq_set, mane_csq), mane_string

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
                csq_string, mane_string = self.get_csq_details(variant)
                candidate_dictionaries.setdefault(variant.sample, []).append(
                    {
                        'variant': self.make_seqr_link(
                            var_string=var_string, sample=sample
                        ),
                        'classes': ', '.join(
                            list(
                                map(
                                    lambda x: self.colours[str(x)],
                                    variant_class_ints,
                                )
                            )
                        ),
                        'symbol': self.panelapp.get(variant.gene).get('symbol'),
                        'csq': csq_string,
                        'mane_select': mane_string,
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
