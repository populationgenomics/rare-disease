"""
A method for constructing a dictionary of CPG sample ID to Seqr Family ID

use the developer console in a browser
load a page in seqr with the network tab open
note the `details` API call
save the response as JSON, and parse the file
pass the result back as a dictionary

Seqr doesn't directly expose an API to access this data

Plan A is contact Seqr directly, so these methods would still be useful
Plan B would be to read similar content from the Sample Metadata DB
Plan ... whatever, would be to do this
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

    # map CPG samples to family ID (unique per project)
    lookup_dict = {}
    for guid, data in details_dict['individualsByGuid'].items():
        if guid in individual_to_cpg.keys():
            lookup_dict[individual_to_cpg[guid]] = data['familyGuid']

    return lookup_dict
