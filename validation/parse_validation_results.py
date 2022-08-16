#!/usr/bin/env python3

"""
run locally to upload the validation results to metamist
"""

import logging
from argparse import ArgumentParser
from csv import DictReader

from cloudpathlib import CloudPath

from sample_metadata.apis import AnalysisApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.model.analysis_status import AnalysisStatus
from sample_metadata.model.analysis_query_model import AnalysisQueryModel


ANAL_API = AnalysisApi()
SUMMARY_KEYS = {
    'TRUTH.TOTAL': 'true_variants',
    'QUERY.TOTAL': 'pipeline_variants',
    'TRUTH.TP': 'matched_variants',
    'TRUTH.FN': 'false_negatives',
    'QUERY.FP': 'false_positives',
    'METRIC.Recall': 'recall',
    'METRIC.Precision': 'precision',
    'METRIC.F1_Score': 'f1_score',
}
VCF_FOLDER = 'single_sample_vcfs'
COMPARISON_FOLDER = 'comparison'

logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


def get_sample_truth(cpg_id: str, dataset: str) -> tuple[str, str]:
    """
    query metamist for the sample truth

    Parameters
    ----------
    cpg_id :
    dataset :

    Returns
    -------
    Path to the truth VCF and corresponding BED file

    """
    a_query_model = AnalysisQueryModel(
        projects=[dataset],
        sample_ids=[cpg_id],
        type=AnalysisType('custom'),
        meta={'type': 'validation'},
        active=True,
    )
    analyses = ANAL_API.query_analyses(analysis_query_model=a_query_model)
    if len(analyses) > 1:
        raise Exception(
            f'Multiple [custom] analysis objects were found for '
            f'{cpg_id}, please set old analyses to active=False'
        )

    if len(analyses) == 0:
        raise Exception(f'{cpg_id} has no Analysis entries')

    analysis_object = analyses[0]
    truth_vcf = analysis_object['output']
    truth_bed = analysis_object['meta'].get('confident_region')
    assert (
        truth_bed and truth_vcf
    ), f'Missing one or both of the truth files: BED: {truth_bed}, VCF: {truth_vcf}'

    return truth_bed, truth_vcf


def post_results(
    cpg_id: str,
    results_path: CloudPath,
    single_sample_vcf: str,
    truth_vcf: str,
    truth_bed: str,
    dataset: str,
):
    """
    parses for all analysis results
    parses the analysis summary
    creates an inactive, completed, QC summary
    posts to metamist

    Parameters
    ----------
    cpg_id :
    results_path :
    single_sample_vcf :
    truth_vcf :
    truth_bed :
    dataset :

    Returns
    -------

    """

    # cast as a list so we can iterate without emptying
    sample_results = list(CloudPath(results_path).glob(f'{cpg_id}*'))

    # pick out the summary file
    summary_file = [file for file in sample_results if 'summary.csv' in file.name][0]

    # populate a dictionary of results for this sample
    summary_data = {
        'type': 'validation',
        'truth_vcf': truth_vcf,
        'truth_bed': truth_bed,
    }
    with summary_file.open() as handle:
        summary_reader = DictReader(handle)
        for line in summary_reader:
            summary_key = f'{line["Type"]}_{line["Filter"]}'
            for sub_key, sub_value in SUMMARY_KEYS.items():
                summary_data[f'{summary_key}::{sub_value}'] = str(line[sub_key])

    for file in sample_results:
        summary_data[file.name.replace(f'{cpg_id}.', '')] = str(file.absolute())

    AnalysisApi().create_new_analysis(
        project=dataset,
        analysis_model=AnalysisModel(
            sample_ids=[cpg_id],
            type=AnalysisType('qc'),
            status=AnalysisStatus('completed'),
            output=single_sample_vcf,
            meta=summary_data,
            active=False,
        ),
    )


def samples_from_vcfs(folder: CloudPath) -> dict[str, str]:
    """

    Parameters
    ----------
    folder :

    Returns
    -------

    """

    vcf_folder = folder / VCF_FOLDER
    sample_vcf_dict = {}
    for file in [file for file in vcf_folder.glob('*') if file.suffix != '.tbi']:
        sample_id = file.name.split('.')[0]
        sample_vcf_dict[sample_id] = str(file.absolute())
    print(sample_vcf_dict)
    return sample_vcf_dict


def main(validation_folder: str, dataset: str):
    """

    Parameters
    ----------
    validation_folder : path to the root of _this_ validation run
    dataset : the dataset to use in metamist

    Returns
    -------

    """

    folder_root = CloudPath(validation_folder)
    sample_dict = samples_from_vcfs(folder_root)
    for cpg_id, ss_vcf in sample_dict.items():
        truth, bed = get_sample_truth(cpg_id=cpg_id, dataset=dataset)
        post_results(
            cpg_id=cpg_id,
            results_path=folder_root / COMPARISON_FOLDER,
            single_sample_vcf=ss_vcf,
            truth_vcf=truth,
            truth_bed=bed,
            dataset=dataset,
        )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        '-i', help='path to the validation report folder (ending in a date)'
    )
    parser.add_argument('-d', help='dataset to use')
    args = parser.parse_args()
    main(args.i, args.d)
