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


DEFAULT_IMAGE = get_config()['workflow']['driver_image']
HAPPY_IMAGE = image_path('happy-vcfeval')
assert DEFAULT_IMAGE, HAPPY_IMAGE

MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')
QUERY_COMPARISON = os.path.join(os.path.dirname(__file__), 'hail_query_validate.py')
OUTPUT_VCFS = output_path('single_sample_vcfs')

# create a logger
logger = logging.getLogger(__file__)


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

    samples_in_jc = set(full_mt.s.collect()).intersection(samples)
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


def comparison_job(batch, ss_vcf: str, sample: str, truth_vcf: str, truth_bed: str):
    """

    Parameters
    ----------
    batch :
    ss_vcf :
    sample :
    truth_vcf :
    truth_bed :

    Returns
    -------

    """

    job = batch.new_job(name=f'Compare {sample}')
    job.image(HAPPY_IMAGE)
    job.memory('20Gi')
    vcf_input = batch.read_input_group(**{'vcf': ss_vcf, 'index': ss_vcf + '.tbi'})
    truth_input = batch.read_input_group(
        **{'vcf': truth_vcf, 'index': truth_vcf + '.tbi'}
    )
    truth_bed = batch.read_input(truth_bed)
    refgenome = batch.read_input(
        'gs://cpg-reference/hg38/v0/dragen_reference/'
        'Homo_sapiens_assembly38_masked.fasta'
    )

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

    # in future don't regenerate SDF...
    job.command(
        f'java -jar -Xmx16G /vcfeval/RTG.jar format -o refgenome_sdf {refgenome} && '
        f'mv {vcf_input["vcf"]} input.vcf.gz && '
        f'mv {vcf_input["index"]} input.vcf.gz.tbi && '
        f'mkdir {job.output} && '
        f'hap.py {truth_input["vcf"]} input.vcf.gz '
        f'-r {refgenome} -R {truth_bed} '
        f'-o {job.output}/output --engine=vcfeval '
        f'--engine-vcfeval-path=/vcfeval/rtg '
        f'--engine-vcfeval-template refgenome.sdf '
        f'--preprocess-truth'
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
    return dict(results)


def get_sample_truth(cpg_id: str) -> tuple[str, str]:
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
        active=True,
    )
    anal_api = AnalysisApi()
    analyses = anal_api.query_analyses(analysis_query_model=a_query_model)
    if len(analyses) != 1:
        logger.error(
            f'Multiple [custom] analysis objects were found for '
            f'{cpg_id}, please set old analyses to active=False'
        )

    analysis_meta = analyses[0]['meta']
    truth_bed = analysis_meta.get('truth_bed')
    truth_vcf = analysis_meta.get('truth_vcf')
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
        name='run validation bits and pieces',
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

    single_sample_files = AnyPath(OUTPUT_VCFS).glob('*.vcf.bgz')
    print(f'Single Sample files: {", ".join(map(str, single_sample_files))}')

    # for each sample, use metamist to pull the corresponding truth and VCF
    # THEN GO AT IT BABY
    # skip any samples without registered truth, complain
    for ss_file in single_sample_files:
        sample_id = ss_file.name.split('.vcf.bgz')[0]
        cpg_id = validation_lookup[sample_id]
        full_path = ss_file.absolute()
        truth_bed, truth_vcf = get_sample_truth(cpg_id)
        logger.info(truth_bed, truth_vcf, full_path, cpg_id, sample_id)
        comparison_job(
            batch=batch,
            ss_vcf=full_path,
            sample=sample_id,
            truth_bed=truth_bed,
            truth_vcf=truth_vcf,
        )

    # twist_bed = 'gs://cpg-validation-test/Twist_Exome_Core_Covered_Targets_hg38.bed'

    batch.run(wait=False)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='input_path')
    parser.add_argument('--header', help='header_lines_file', default=None)
    args = parser.parse_args()
    main(input_file=args.i, header=args.header)
