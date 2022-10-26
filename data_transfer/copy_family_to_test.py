#!/usr/bin/env python3

"""

Given a project ID and one or more family IDs, copies cram and gvcf
files for each individual in the family into the project test bucket.

"""

from argparse import ArgumentParser
from collections import defaultdict
import logging
import sys
import subprocess


from sample_metadata.apis import FamilyApi, AnalysisApi, SampleApi
from sample_metadata.models import BodyGetSamples, AnalysisType


def get_family_id_to_participant_map(project: str) -> dict[str, list[dict]]:
    """
    Generate a map of external family IDs to a list of internal participant IDs
    """

    result = defaultdict(list)
    pedigree = FamilyApi().get_pedigree(
        project=project,
        replace_with_participant_external_ids=False,
        replace_with_family_external_ids=True,
    )

    for individual in pedigree:
        result[individual['family_id']].append(individual)

    return dict(result)


def copy_to_test(project: str, path: str):
    """
    Copy a single file from a main bucket path to the equivalent test bucket
    """
    test_path = path.replace(
        f'cpg-{project}-main',
        f'cpg-{project}-test',
    )

    subprocess.run(f'gsutil cp {path} {test_path}', shell=True, check=True)
    logging.info(f'Copied {test_path}')


def main(
    project: str,
    family_ids: list[str],
):
    """

    Parameters
    ----------
    project : meatamist project name
    family_ids : a list of external family IDs to transfer to test

    Returns
    -------

    """

    # Find all participants in the nominated families
    participants_by_ext_family_id = get_family_id_to_participant_map(project)
    participant_ids = []
    unknown_family_id = False
    for family_id in family_ids:
        try:
            for participant in participants_by_ext_family_id[family_id]:
                participant_ids.append(participant['individual_id'])
        except KeyError:
            unknown_family_id = True
            print(
                f'Error: "{family_id}" is not a valid external family ID in '
                f'project "{project}".',
                file=sys.stderr,
            )

    if unknown_family_id:
        sys.exit(1)

    # Retrieve active samples for these participants
    samples = SampleApi().get_samples(
        body_get_samples=BodyGetSamples(
            project_ids=[project], participant_ids=participant_ids, active=True
        )
    )

    # Retrieve latest crams and gvcfs for selected samples
    sample_ids = [sample['id'] for sample in samples]
    latest_crams = AnalysisApi().get_latest_analysis_for_samples_and_type(
        AnalysisType('cram'), project, request_body=sample_ids
    )
    latest_gvcfs = AnalysisApi().get_latest_analysis_for_samples_and_type(
        AnalysisType('gvcf'), project, request_body=sample_ids
    )

    # Copy files to test
    for cram in latest_crams:
        copy_to_test(project, cram['output'])
    for gvcf in latest_gvcfs:
        copy_to_test(project, gvcf['output'])


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%M-%d %H:%M:%S',
        stream=sys.stderr,
    )

    parser = ArgumentParser()
    parser.add_argument('-p', '--project', help='Project name', required=True)
    parser.add_argument(
        'family_ids', nargs='+', help='External family IDs to be transferred.'
    )

    args, unknown = parser.parse_known_args()

    if unknown:
        raise Exception(f'Unknown args, could not parse: "{unknown}"')

    # init_batch()

    main(
        project=args.project,
        family_ids=args.family_ids,
    )
