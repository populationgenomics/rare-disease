#!/usr/bin/env python3

"""
run locally to upload the validation results to metamist
"""

from argparse import ArgumentParser
from csv import DictReader

from cloudpathlib import CloudPath

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path

from sample_metadata.apis import AnalysisApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.model.analysis_status import AnalysisStatus


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


def main(cpg_id: str, single_sample_vcf: str, truth_vcf: str, truth_bed: str):
    """
    parses for this sample's analysis results
    parses the analysis summary CSV
    posts a QC entry to metamist with meta.type=validation_result

    Parameters
    ----------
    cpg_id : CPG#### Identifier
    single_sample_vcf : The specific VCF generated and used in validation
    truth_vcf : the sample truth VCF
    truth_bed : the confident regions BED
    """

    # cast as a list so we can iterate without emptying
    sample_results = list(CloudPath(output_path('comparison')).glob(f'{cpg_id}*'))

    # pick out the summary file
    summary_file = [file for file in sample_results if 'summary.csv' in file.name][0]

    # populate a dictionary of results for this sample
    summary_data = {
        'type': 'validation_result',
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

    # should the output be the summary file rather than the VCF? Hmm
    AnalysisApi().create_new_analysis(
        project=get_config()['workflow']['dataset'],
        analysis_model=AnalysisModel(
            sample_ids=[cpg_id],
            type=AnalysisType('qc'),
            status=AnalysisStatus('completed'),
            output=single_sample_vcf,
            meta=summary_data,
            active=True,
        ),
    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--id', help='CPG ID for this sample')
    parser.add_argument('--ss', help='single sample VCF used')
    parser.add_argument('-t', help='truth VCF used')
    parser.add_argument('-b', help='BED used')
    args = parser.parse_args()
    main(args.id, single_sample_vcf=args.ss, truth_vcf=args.t, truth_bed=args.b)
