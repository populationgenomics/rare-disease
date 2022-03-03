"""
generate a ped file on the fly using the sample-metadata api client
additional (optional) argument will remove all family associations
"""


from typing import Dict, List
from sample_metadata.apis import FamilyApi, ParticipantApi

import click


PED_KEYS = [
    'family_id',
    'individual_id',
    'paternal_id',
    'maternal_id',
    'sex',
    'affected',
]


def get_clean_pedigree(
    pedigree_dicts: List[Dict[str, str]],
    sample_to_cpg_dict: Dict[str, str],
    singles: bool,
) -> List[Dict[str, str]]:
    """
    maybe use dictwriter?
    """

    new_entries = []

    for ped_entry in pedigree_dicts:

        if ped_entry['individual_id'] not in sample_to_cpg_dict:
            continue

        # update the sample IDs
        ped_entry['individual_id'] = sample_to_cpg_dict.get(ped_entry['individual_id'])

        # remove parents and assign an individual sample ID
        if singles:
            ped_entry['paternal_id'] = ''
            ped_entry['maternal_id'] = ''
            ped_entry['family_id'] = ped_entry['individual_id']

        else:
            ped_entry['paternal_id'] = sample_to_cpg_dict.get(
                ped_entry['paternal_id'], ''
            )
            ped_entry['maternal_id'] = sample_to_cpg_dict.get(
                ped_entry['maternal_id'], ''
            )

        new_entries.append(ped_entry)

    return new_entries


def write_pedigree(clean_pedigree: List[Dict[str, str]], output: str):
    """
    stick it in a file
    """
    with open(output, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(PED_KEYS) + '\n')
        for entry in clean_pedigree:
            handle.write('\t'.join([str(entry.get(key)) for key in PED_KEYS]) + '\n')


def get_pedigree_for_project(project: str) -> List[Dict[str, str]]:
    """
    fetches the project pedigree from sample-metadata
    one list entry for each participant
    :param project:
    :return:
    """
    return FamilyApi().get_pedigree(project=project)


def ext_to_int_sample_map(project: str) -> Dict[str, str]:
    """
    fetches the mapping dictionary, so external IDs can be
    translated to the corresponding CPG ID
    :param project:
    :return:
    """
    return ParticipantApi().get_external_participant_id_to_internal_sample_id(
        project=project
    )


@click.command()
@click.option(
    '--project',
    help='the name of the project to use in API queries',
)
@click.option(
    '--singles',
    default=False,
    is_flag=True,
    help='remake the pedigree as singletons',
)
@click.option(
    '--output',
    help='write the new PED file here',
)
def main(project: str, singles: bool, output: str):
    """

    :param project: may be able to retrieve this from the environment
    :param singles: whether to split the pedigree(s) into singletons
    :param output: path to write new PED file to
    """

    # get the list of all pedigree members as list of dictionaries
    pedigree_dicts = get_pedigree_for_project(project=project)

    # endpoint gives list of lists e.g. [['A1234567_proband', 'CPG12341']]
    sample_to_cpg_dict = ext_to_int_sample_map(project=project)

    clean_pedigree = get_clean_pedigree(
        pedigree_dicts=pedigree_dicts,
        sample_to_cpg_dict=sample_to_cpg_dict,
        singles=singles,
    )

    write_pedigree(clean_pedigree, output)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
