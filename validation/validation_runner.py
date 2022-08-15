#!/usr/bin/env python3

"""
wraps the validation script(s)
"""

import logging
import os
from pathlib import Path
from argparse import ArgumentParser

from cloudpathlib import AnyPath
import hail as hl
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    init_batch,
    output_path,
    remote_tmpdir,
    image_path,
)
from sample_metadata.apis import AnalysisApi, ParticipantApi
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_type import AnalysisType


MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')
OUTPUT_VCFS = output_path('single_sample_vcfs')
ANAL_API = AnalysisApi()

# create a logger
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


def mt_to_vcf(input_mt: str, header_lines: str | None, samples: set[str]):
    """
    takes a MT and converts to VCF
    adds in extra header lines for VQSR filters
    :param input_mt:
    :param header_lines:
    :param samples:
    :return:
    """
    init_batch()

    full_mt = hl.read_matrix_table(input_mt)
    all_jc_samples = full_mt.s.collect()
    logging.info(f'All Joint-call samples: {all_jc_samples}')

    samples_in_jc = set(all_jc_samples).intersection(samples)
    logging.info(
        f'Extracting these samples from the joint-call: {" ".join(samples_in_jc)}'
    )

    # extract all present samples into a separate file
    for sample in samples_in_jc:

        logging.info(f'Processing: {sample}')

        sample_path = os.path.join(OUTPUT_VCFS, f'{sample}.vcf.bgz')

        if AnyPath(sample_path).exists():
            print(f'no action taken, {sample_path} already exists')
            continue

        # filter to this column
        mt = full_mt.filter_cols(full_mt.s == sample)

        # filter to this sample's non-ref calls
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)
        hl.export_vcf(
            mt,
            sample_path,
            append_to_header=header_lines,
            tabix=True,
        )


def comparison_job(
    batch, ss_vcf: str, sample: str, truth_vcf: str, truth_bed: str, reference_sdf: str
):
    """

    Parameters
    ----------
    batch :
    ss_vcf :
    sample :
    truth_vcf :
    truth_bed :
    reference_sdf :

    Returns
    -------

    """

    job = batch.new_job(name=f'Compare {sample}')
    job.image(image_path('happy'))
    job.memory('20Gi')
    job.storage('20Gi')
    vcf_input = batch.read_input_group(**{'vcf': ss_vcf, 'index': f'{ss_vcf}.tbi'})
    truth_input = batch.read_input_group(
        **{'vcf': truth_vcf, 'index': f'{truth_vcf}.tbi'}
    )
    truth_bed = batch.read_input(truth_bed)
    refgenome = (
        'gs://cpg-reference/hg38/v0/dragen_reference/'
        'Homo_sapiens_assembly38_masked.fasta'
    )
    batch_ref = batch.read_input_group(
        **{'fasta': refgenome, 'index': f'{refgenome}.fai'}
    )

    print(reference_sdf)
    # # sdf loading as a Glob operation
    # sdf = batch.read_input_group(
    #     **{file.name: file.as_uri() for file in AnyPath(reference_sdf).glob('*')}
    # )

    # hap.py outputs:
    # output.extended.csv
    # output.roc.all.csv.gz
    # output.metrics.json.gz
    # output.runinfo.json
    # output.roc.Locations.INDEL.PASS.csv.gz
    # output.summary.csv
    # output.roc.Locations.INDEL.csv.gz
    # output.vcf.gz
    # output.roc.Locations.SNP.PASS.csv.gz
    # output.vcf.gz.tbi
    # output.roc.Locations.SNP.csv.gz

    job.declare_resource_group(
        output={
            'happy_extended.csv': '{root}/output.extended.csv',
            'happy.vcf.gz': '{root}/output.vcf.gz',
            'happy.vcf.gz.tbi': '{root}/output.vcf.gz.tbi',
            'happy_roc.all.csv.gz': '{root}/output.roc.all.csv.gz',
            'happy_metrics.json.gz': '{root}/output.metrics.json.gz',
            'happy_runinfo.json': '{root}/output.runinfo.json',
            'summary.csv': '{root}/output.summary.csv',
        }
    )

    # pre.py use is cancelled out for now
    # # set arguments for pre-py
    pre_args = f'-r {batch_ref["fasta"]} --pass-only -R {truth_bed}'

    # pre-process the truth and test data
    job.command(
        f'pre.py {pre_args} {truth_input["vcf"]} truth.vcf.gz && '
        f'pre.py {pre_args} {vcf_input["vcf"]} query.vcf.gz && '
        f'mkdir {job.output} && '
        f'hap.py truth.vcf.gz query.vcf.gz '
        f'-r {batch_ref["fasta"]} -R {truth_bed} '
        f'-o {job.output}/output --threads 10 '
    )
    batch.write_output(job.output, os.path.join(output_path('comparison'), sample))


def get_validation_samples() -> dict[str, str]:
    """
    query metamist for all sample IDs
    return a dict of {external: CPG ID}
    Returns
    -------

    """

    party = ParticipantApi()
    results = party.get_external_participant_id_to_internal_sample_id(
        project=get_config()['workflow']['dataset']
    )
    # pylint: disable=unnecessary-comprehension
    return {y: x for x, y in results}


def get_sample_truth(cpg_id: str) -> tuple[str | None, str | None]:
    """
    query metamist for the sample truth
    Parameters
    ----------
    cpg_id :

    Returns
    -------
    Path to the truth VCF and corresponding BED file

    """
    a_query_model = AnalysisQueryModel(
        projects=[get_config()['workflow']['dataset']],
        sample_ids=[cpg_id],
        type=AnalysisType('custom'),
        meta={'type': 'validation'},
        active=True,
    )
    analyses = ANAL_API.query_analyses(analysis_query_model=a_query_model)
    if len(analyses) > 1:
        logger.error(
            f'Multiple [custom] analysis objects were found for '
            f'{cpg_id}, please set old analyses to active=False'
        )
        return None, None

    if len(analyses) == 0:
        logger.error(f'{cpg_id} has no Analysis entries')
        return None, None

    analysis_object = analyses[0]
    truth_vcf = analysis_object['output']
    truth_bed = analysis_object['meta'].get('confident_region')
    assert (
        truth_bed and truth_vcf
    ), f'Missing one or both of the truth files: BED: {truth_bed}, VCF: {truth_vcf}'

    return truth_bed, truth_vcf


def main(input_file: str, header: str | None):
    """

    Parameters
    ----------
    input_file :
    header :

    Returns
    -------

    """

    validation_lookup = get_validation_samples()

    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='run without vcfeval preprocessing',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory='highmem',
    )

    input_path = Path(input_file)
    if input_path.suffix != '.mt':
        logger.error('Expected a MT as input file')
        raise Exception('Expected a MT as input file')

    # requires conversion to vcf
    mt_to_vcf(
        input_mt=input_file,
        header_lines=header,
        samples=set(validation_lookup.keys()),
    )

    single_sample_files = list(AnyPath(OUTPUT_VCFS).glob('*.vcf.bgz'))

    ref_sdf = 'gs://cpg-validation-test/refgenome_sdf'

    # for each sample, use metamist to pull the corresponding truth and VCF
    # THEN GO AT IT BABY
    # skip any samples without registered truth, complain
    scheduled_jobs = False
    for ss_file in single_sample_files:
        cpg_id = ss_file.name.split('.vcf.bgz')[0]
        full_path = ss_file.absolute()
        truth_bed, truth_vcf = get_sample_truth(cpg_id)
        if truth_bed is None:
            logger.error(f'Skipping validation run for {cpg_id}')
            continue
        comparison_job(
            batch=batch,
            ss_vcf=str(full_path),
            sample=cpg_id,
            truth_bed=str(truth_bed),
            truth_vcf=str(truth_vcf),
            reference_sdf=ref_sdf,
        )
        scheduled_jobs = True

    # twist_bed = 'gs://cpg-validation-test/Twist_Exome_Core_Covered_Targets_hg38.bed'

    if scheduled_jobs:
        batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('-i', help='input_path')
    parser.add_argument('--header', help='header_lines_file', default=None)
    args = parser.parse_args()
    main(input_file=args.i, header=args.header)
