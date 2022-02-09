#!usr/bin/env python3

"""
little script, takes three inputs
- PED
- PanelApp data
- Post-Slivar VCF

Reads the VCF, grouping variants by gene
Reads the inheritance pattern to use from panelapp
Reports variants which fit the given pattern
Uses the PED file to:
    - appropriately classify X-chrom inheritance
    - avoid reporting on any unaffected participants

Crucially, the PED is only used in the context of singleton
analysis here, so we are not checking family members
"""
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional
from itertools import combinations
from csv import DictReader
import json
import logging
import os
import click

from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant


CHROMOSOMES = [f"chr{x}" for x in list(range(1, 23)) + ['X', 'Y', "M"]]


@dataclass
class PedPerson:
    """
    holds attributes about a single PED file entry
    """

    sample: str
    male: bool
    affected: bool


def parse_ped_simple(ped: str) -> Dict[str, PedPerson]:
    """
    just take individual attributes - sample ID, sex, affected
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


def variants_from_chunk(vcf: str, chrom: str):
    """
    takes the requested chromosome, and returns all variants by gene
    :param vcf: path to the vcf
    :param chrom: the CHR to limit this search to
    :return:
    """

    vcf_dict = {}
    reader = VCF(vcf, strict_gt=True)
    for var in reader(chrom):
        vcf_dict.setdefault(var.INFO.get('vep_gene_ids'), []).append(var)
    return vcf_dict


@click.command()
@click.option("--vcf", "vcf", help="VCF to parse")
@click.option("--panel", "panel", help="JSON file containing the PanelApp MOI digest")
@click.option(
    "--ped", "ped", help="PED file, not really required for singleton", default=None
)
def main(vcf: str, panel: str, ped: Optional[str] = None):
    """

    :param vcf:
    :param panel:
    :param ped: used to limit inheritance to affected
    :return:
    """
    assert os.path.exists(panel), f'PanelData required, {panel} doesn\'t exist'
    assert os.path.exists(vcf), f'VCF required, {vcf} doesn\'t exist'
    if ped is not None:
        assert os.path.exists(ped), f'PED required if specified, {ped} doesn\'t exist'
        ped = parse_ped_simple(ped)

    with open(panel, 'r', encoding='utf-8') as handle:
        panel_dict = json.load(handle)

    valid_variants = {}
    for chrom in CHROMOSOMES:
        chrom_vars = variants_from_chunk(vcf, chrom)
        for ensg, varlist in chrom_vars.items():

            for key, value in check_moi(
                varlist=varlist, inheritance=panel_dict.get(ensg), ensg=ensg, ped=ped
            ).items():
                valid_variants.setdefault(key, []).extend(value)

    for sample, patterns in valid_variants.items():

        # if we had a pedigree, don't print out data on
        # unaffected samples (parents in this joint call)
        if ped is not None and not ped.get(sample).affected:
            continue
        print(f'{sample}:')
        for pat in patterns:
            print(f'\t{pat}\n')


def check_moi(
    varlist: List[Variant],
    inheritance: str,
    ensg: str,
    ped: Optional[Dict[str, PedPerson]] = None,
) -> Dict[str, List[str]]:
    """
    simplify all this crap, probably a helper method
    :param varlist: list of variants within this gene
    :param inheritance: the simplified pattern from PanelApp
    :param ensg:
    :param ped:
    :return:
    """
    good_fits = {}

    if inheritance in {'Monoallelic', 'Mono_And_Biallelic', 'Hemi_Mono_In_Female'}:
        for var in varlist:
            for sample in [
                sam for sam in var.INFO.get('het', '').split(',') if sam != ''
            ]:
                good_fits.setdefault(sample, []).append(
                    f'{string_repr(var)},0/1,{ensg}:{inheritance}'
                )
            for sample in [
                sam for sam in var.INFO.get('hom', '').split(',') if sam != ''
            ]:
                good_fits.setdefault(sample, []).append(
                    f'{string_repr(var)},1/1,{ensg}:{inheritance}'
                )

    elif inheritance in {'Biallelic', 'Hemi_Bi_In_Female'}:
        for var in varlist:

            # hemi-het is valid for males
            for sample in [
                sam
                for sam in var.INFO.get('het', '').split(',')
                if sam != ''
                and ped is not None
                and ped.get(sam).male
                and inheritance == 'Hemi_Bi_In_Female'
            ]:
                good_fits.setdefault(sample, []).append(
                    f'{string_repr(var)},0/1,{ensg}:{inheritance}::Male'
                )
            for sample in [
                sam for sam in var.INFO.get('hom', '').split(',') if sam != ''
            ]:
                good_fits.setdefault(sample, []).append(
                    f'{string_repr(var)},1/1,{ensg}:{inheritance}'
                )
        # iterate over all variant combinations
        # each pair considered only once, orientation not important
        for var1, var2 in combinations(varlist, 2):

            # if either one of these variants lacks HET calls, goto next pair
            if var1.INFO.get('het') is None or var2.INFO.get('het') is None:
                continue

            # for each sample in the intersection of these two het groups
            for sample in set(var1.INFO.get('het').split(',')).intersection(
                set(var2.INFO.get('het').split(','))
            ):
                good_fits.setdefault(sample, []).append(
                    f'{string_repr(var1)},{string_repr(var2)},{ensg}:{inheritance}\tHET'
                )

    return good_fits


@lru_cache()
def string_repr(var: Variant):
    """
    just a string printer
    :param var:
    :return:
    """
    return f'{var.CHROM}:{var.POS}:{var.REF}:{var.ALT[0]}'


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)
    main()  # pylint: disable=E1120
