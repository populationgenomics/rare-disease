"""
simplify the results of the prev. script into a TSV
"""

import json

INPUT = 'this_local_filename.json'
OUTPUT = 'simplified.tsv'

with open(INPUT) as handle:
    data = json.load(handle)

gene_dict: dict = {}

for content in data.values():
    for gene_name, gene_data in content['genes'].items():
        if gene_data['confidence_level'] != '3':
            continue

        gene_dict.setdefault(
            gene_name,
            {'pmids': set(), 'disorders': set(), 'ensg': '', 'phenotypes': set()},
        )
        gene_dict[gene_name]['phenotypes'].update(gene_data['phenotypes'])
        gene_dict[gene_name]['pmids'].update(gene_data['publications'])
        gene_dict[gene_name]['disorders'].update(content['disorders'])
        gene_dict[gene_name]['ensg'] = gene_data['ensembl_id']

with open(OUTPUT, 'w', encoding='utf-8') as handle:
    handle.write('gene\tpmids\tdisorders\tphenotypes\tensg\n')
    for gene_name, gene_data in gene_dict.items():
        handle.write(
            f'{gene_name}\t'
            f'{",".join(gene_data["pmids"])}\t'
            f'{",".join(gene_data["disorders"])}\t'
            f'{",".join(gene_data["phenotypes"])}\t'
            f'{gene_data["ensg"]}\n',
        )
