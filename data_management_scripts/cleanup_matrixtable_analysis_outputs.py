#!/usr/bin/env python3  # noqa: EXE001
"""
  !!    THIS SCRIPT CAN DELETE PRODUCTION DATA. USE WITH CAUTION.     !!

Finds Metamist analyses with matrixtable outputs and deletes all but the most recent one.
Separates the analyses by dataset, sequencing type, and analysis type, keeping the most recent one for each.

This script is intended to be run regularly to clean up old matrixtables based on Metamist completion timestamps.
"""

import logging
import sys
from datetime import datetime

import click
from cloudpathlib import CloudPath
from metamist.graphql import gql, query

MT_ANALYSIS_STAGE_MAP = {
    'custom': 'AnnotateDataset',
    'sv': 'AnnotateDatasetSv',
}

MT_ANALYSES_QUERY = gql(
    """
    query DatasetMtAnalyses($dataset: String!, $analysisType: String!) {
        project(name: $dataset) {
            analyses(type: {eq: $analysisType}) {
                id
                output
                timestampCompleted
                meta
            }
        }
    }
    """,
)


def get_mt_analyses_for_dataset(
    dataset: str,
    sequencing_type: str,
    mt_analysis_type: str,
):
    """Finds the 'custom' OR 'sv' matrix table analyses for a dataset and sequencing type"""
    logging.getLogger().setLevel(logging.WARN)
    dataset_mt_analyses = query(
        MT_ANALYSES_QUERY,
        variables={'dataset': dataset, 'analysisType': mt_analysis_type},
    )
    logging.getLogger().setLevel(logging.INFO)

    analyses = []
    for analysis in dataset_mt_analyses['project']['analyses']:
        if (
            analysis['meta'].get('stage') == MT_ANALYSIS_STAGE_MAP[mt_analysis_type]
            and analysis['meta'].get('sequencing_type') == sequencing_type
            and analysis['output'].endswith('.mt')
        ):
            analyses.append(analysis)

    return analyses


def analyses_to_delete(analyses: list[dict[str, str]]):
    """Takes the list of analyses and keeps only the most recently completed one, returning the rest"""
    analyses = sorted(
        analyses,
        key=lambda x: datetime.fromisoformat(x['timestampCompleted']),
    )

    logging.info(
        f'Found latest analysis: \t ID: {analyses[-1]["id"]} \t Output: {analyses[-1]["output"]} \t Timestamp: {analyses[-1]["timestampCompleted"]}',
    )
    logging.info(f'Found {len(analyses[:-1])} older analyses:')
    for analysis in reversed(analyses[:-1]):
        logging.info(
            f'ID: {analysis["id"]} \t Output: {analysis["output"]} \t Timestamp: {analysis["timestampCompleted"]}',
        )

    return analyses[:-1]


@click.command()
@click.option(
    '--dry-run',
    is_flag=True,
    help='Prints the list of matrix tables to delete without deleting them',
)
@click.option(
    '--datasets',
    '-d',
    multiple=True,
    help='The datasets to check for matrix tables',
)
@click.option(
    '--sequencing-types',
    '-s',
    multiple=True,
    help='The sequencing types to check for matrix tables',
    default=['genome', 'exome'],
)
@click.option(
    '--mt-analysis-types',
    '-t',
    multiple=True,
    help='The analysis types to check for matrix tables',
    default=list(MT_ANALYSIS_STAGE_MAP.keys()),
)
def main(
    dry_run: bool,
    datasets: list[str],
    sequencing_types: list[str],
    mt_analysis_types: list[str],
):
    """
    Finds all Metamist analysis entries corresponding to dataset .mt files for the input datasets,
    sequencing types, and analysis types - 'custom' or 'sv' (as of 17/11/2023)
    Deletes all but the most recent matrix table for each dataset, sequencing type, and analysis type.
    """
    if any(
        mt_analysis_type not in MT_ANALYSIS_STAGE_MAP.keys()
        for mt_analysis_type in mt_analysis_types
    ):
        raise ValueError(f'Invalid mt_analysis_types: {mt_analysis_types}')

    matrixtables_to_delete = []
    for dataset in datasets:
        logging.info(f'Checking {dataset}')
        for sequencing_type in [st.lower() for st in sequencing_types]:
            for mt_analysis_type in [t.lower() for t in mt_analysis_types]:
                logging.info(
                    f'{dataset} :: {sequencing_type.upper()} :: Checking for analysis of type {mt_analysis_type}',
                )

                analyses = get_mt_analyses_for_dataset(
                    dataset,
                    sequencing_type,
                    mt_analysis_type,
                )
                logging.info(
                    f'{dataset} :: {sequencing_type.upper()} :: Found {len(analyses)} {mt_analysis_type} analyses',
                )
                if not analyses:
                    continue
                analyses_to_remove = analyses_to_delete(analyses)
                logging.info(
                    f'{dataset}:: {sequencing_type.upper()} :: {mt_analysis_type} analysis outputs to delete:',
                )

                matrixtables_to_delete.extend(analyses_to_remove)

    if dry_run:
        logging.info('<< Dry run, not deleting anything >>')
        logging.info(
            f'Would have deleted {len(matrixtables_to_delete)} .mt analysis outputs from {len(datasets)} datasets',
        )

    else:
        for analysis in matrixtables_to_delete:
            mt = CloudPath(analysis['output'])
            logging.info(f'Deleting {analysis["output"]}')
            mt.rmtree()  # use rmtree to delete the .mt folder and all its contents

        logging.info(
            f'Deleted {len(matrixtables_to_delete)} .mt analysis outputs from {len(datasets)} datasets',
        )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=no-value-for-parameter
