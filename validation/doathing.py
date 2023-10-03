#!/usr/bin/env python3


"""
Handles the validation process
Input:
- Path to MatrixTable
Process:
- Query Metamist for all validation samples
- Extract each validation sample in the MT into a separate VCF
- For each extracted VCF, find matched truth from Metamist
- Run Hap.py validation, using vcfeval normalisation engine
- Log results of validation into Metamist
Results:
- Single Sample VCF
- Per-sample entry in Metamist logging results
"""


import logging
import os
from argparse import ArgumentParser
from pathlib import Path

from cpg_workflows.batch import get_batch
from sample_metadata.apis import AnalysisApi, ParticipantApi
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_type import AnalysisType

MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'parse_validation_results.py')
REF_SDF = 'gs://cpg-validation-test/refgenome_sdf'

# create a logger
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


def get_validation_samples() -> dict[str, str]:
    """
    query metamist for all sample IDs in this dataset

    Returns
    -------
    dict of {external: CPG ID}
    """

    party = ParticipantApi()
    results = party.get_external_participant_id_to_internal_sample_id(
        project='validation',
    )
    return {y: x for x, y in results}


def get_sample_truth(cpg_id: str) -> tuple[str, str]:
    """
    query metamist for the sample truth

    Parameters
    ----------
    cpg_id : CPG### format ID

    Returns
    -------
    Path to the truth VCF and confident region BED file
    """

    a_query_model = AnalysisQueryModel(
        projects=['validation'],
        sample_ids=[cpg_id],
        type=AnalysisType('custom'),
        meta={'type': 'validation'},
        active=True,
    )
    analyses = AnalysisApi().query_analyses(analysis_query_model=a_query_model)
    print(analyses)
    if len(analyses) > 1:
        raise Exception(
            f'Multiple [custom] analysis objects were found for '
            f'{cpg_id}, please set old analyses to active=False',
        )

    if len(analyses) == 0:
        raise Exception(f'{cpg_id} has no Analysis entries')

    analysis_object = analyses[0]
    truth_vcf = analysis_object.get('output')
    truth_bed = analysis_object['meta'].get('confident_region')
    assert (
        truth_bed and truth_vcf
    ), f'Missing one or both of the truth files: BED: {truth_bed}, VCF: {truth_vcf}'

    return truth_vcf, truth_bed


def main(input_file: str, stratification: str | None, dry_run: bool = False):
    """
    Parameters
    ----------
    input_file : path to the MT representing this joint-call
    stratification : the path to the stratification BED files
    dry_run : if True, prevent writes to Metamist
    """

    input_path = Path(input_file)
    if input_path.suffix != '.mt':
        logger.error('Expected a MT as input file')
        raise Exception('Expected a MT as input file')

    # # set the path for this output
    # process the MT to get the name

    validation_lookup = get_validation_samples()
    truth_vcf, truth_bed = get_sample_truth('cpg')
