"""
This script is home to all the content required to parse Seqr metadata
This is unlikely to be deployed in anything resembling this format

Plan A is contact Seqr directly, so these methods would still be useful
Plan B would be to read similar content from the Sample Metadata DB
Plan ... whatever, would be to do this

_THIS_ is using the developer console in the browser to monitor the data
loaded by Seqr, save the response data, and parse it in static form
dumping the result out as a dictionary
"""

from typing import Dict
import json
import os


def get_seqr_details(seqr_meta: str) -> Dict[str, Dict[str, str]]:
    """
    this would be substituted for a GET call
    :return: dict of the semi-processed output
    """

    if not os.path.exists(seqr_meta):
        return {}

    with open(seqr_meta, 'r', encoding='utf-8') as handle:
        details_dict = json.load(handle)

    # map CPG ID to individual GUID
    individual_to_cpg = {
        sample['individualGuid']: sample['sampleId']
        for seqr_sample_id, sample in details_dict['samplesByGuid'].items()
    }

    # map CPG samples to families
    # the family ID itself is unique to a project
    # meaning that project doesn't need to feature in the URL
    lookup_dict = {}
    for guid, data in details_dict['individualsByGuid'].items():
        if guid in individual_to_cpg.keys():
            lookup_dict[individual_to_cpg[guid]] = {
                'family': data['familyGuid'],
            }

    return lookup_dict