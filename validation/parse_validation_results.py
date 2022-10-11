#!/usr/bin/env python3

"""
Runs after the validation process to digest results
A summary of those results are posted into Metamist
"""


import logging
from argparse import ArgumentParser
from csv import DictReader
import json

from cpg_utils import to_path
from cpg_utils.config import get_config

from sample_metadata.apis import AnalysisApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.model.analysis_status import AnalysisStatus
from sample_metadata.model.analysis_query_model import AnalysisQueryModel


ANAL_API = AnalysisApi()
SUMMARY_KEYS = {
    'TRUTH.TOTAL': 'true_variants',
    'METRIC.Recall': 'recall',
    'METRIC.Precision': 'precision',
}


def check_for_prior_result(cpg_id: str, comparison_folder: str) -> bool:
    """
    queries Metamist to make sure that the results were not already
    logged in a prior process

    Parameters
    ----------
    cpg_id : str
    comparison_folder : str

    Returns
    -------
    bool - if True, this result set was already logged
    """

    a_query_model = AnalysisQueryModel(
        projects=[get_config()['workflow']['dataset']],
        sample_ids=[cpg_id],
        output=comparison_folder,
    )
    analyses = ANAL_API.query_analyses(analysis_query_model=a_query_model)
    return len(analyses) > 0


def main(
    cpg_id: str,
    comparison_folder: str,
    single_sample_vcf: str,
    truth_vcf: str,
    truth_bed: str,
    joint_mt: str,
    stratified: str | None = None,
    dry_run: bool = False,
):
    """
    parses for this sample's analysis results
    parses the analysis summary CSV
    posts a QC entry to metamist with meta.type=validation_result

    Parameters
    ----------
    cpg_id : CPG#### Identifier
    comparison_folder :
    single_sample_vcf : The specific VCF generated and used in validation
    truth_vcf : the sample truth VCF
    truth_bed : the confident regions BED
    joint_mt : the original joint-call MatrixTable
    stratified : if analysis was stratified, this is the location of files
    dry_run : don't post to metamist
    """

    # if a result already exists, quietly exit so as not to cancel other sample's jobs
    if check_for_prior_result(cpg_id=cpg_id, comparison_folder=comparison_folder):
        logging.info(
            f'Sample {cpg_id} already has validation '
            f'results logged from {comparison_folder}, quitting'
        )
        return

    # list from the generator, so we can re-iterate
    sample_results = list(to_path(comparison_folder).glob(f'{cpg_id}*'))

    # pick out the summary file
    summary_file = [file for file in sample_results if 'extended.csv' in file.name][0]

    # populate a dictionary of results for this sample
    summary_data = {
        'type': 'validation_result',
        'source_mt': joint_mt,
        'query_vcf': single_sample_vcf,
        'truth_vcf': truth_vcf,
        'truth_bed': truth_bed,
    }

    if stratified:
        summary_data['stratified'] = stratified

    with summary_file.open() as handle:
        summary_reader = DictReader(handle)
        for line in summary_reader:
            if line['Filter'] != 'PASS' or line['Subtype'] != '*':
                continue

            summary_key = f'{line["Type"]}_{line["Subset"]}'
            for sub_key, sub_value in SUMMARY_KEYS.items():
                summary_data[f'{summary_key}::{sub_value}'] = str(line[sub_key])

    # store the full paths of all files created during the analysis
    for file in sample_results:
        summary_data[file.name.replace(f'{cpg_id}.', '')] = str(file.absolute())

    # print the contents, even if the metamist write fails (e.g. on test)
    logging.info(json.dumps(summary_data, indent=True))

    if not dry_run:
        AnalysisApi().create_new_analysis(
            project=get_config()['workflow']['dataset'],
            analysis_model=AnalysisModel(
                sample_ids=[cpg_id],
                type=AnalysisType('qc'),
                status=AnalysisStatus('completed'),
                output=comparison_folder,
                meta=summary_data,
                active=True,
            ),
        )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('--id', help='CPG ID for this sample')
    parser.add_argument('--folder', help='Location for results')
    parser.add_argument('--ss', help='single sample VCF used')
    parser.add_argument('-t', help='truth VCF used')
    parser.add_argument('-b', help='BED used')
    parser.add_argument('--mt', help='Multisample MT')
    parser.add_argument(
        '--stratified', help='Stratification files, if used', required=False
    )
    parser.add_argument(
        '--dry_run', action='store_true', help="if present, don't write to metamist"
    )
    args = parser.parse_args()
    main(
        cpg_id=args.id,
        comparison_folder=args.folder,
        single_sample_vcf=args.ss,
        truth_vcf=args.t,
        truth_bed=args.b,
        joint_mt=args.mt,
        stratified=args.stratified,
        dry_run=args.dry_run,
    )
