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
from pathlib import Path
from argparse import ArgumentParser

import hailtop.batch as hb
from cloudpathlib import AnyPath
import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.git import (
    prepare_git_job,
    get_git_commit_ref_of_current_repository,
    get_organisation_name_from_current_directory,
    get_repo_name_from_current_directory,
)
from cpg_utils.hail_batch import (
    init_batch,
    output_path,
    remote_tmpdir,
    image_path,
    copy_common_env,
    authenticate_cloud_credentials_in_job,
)
from cpg_utils.workflows.batch import Batch
from sample_metadata.apis import AnalysisApi, ParticipantApi
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_type import AnalysisType


MT_TO_VCF_SCRIPT = os.path.join(os.path.dirname(__file__), 'mt_to_vcf.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'parse_validation_results.py')
REF_SDF = 'gs://cpg-validation-test/refgenome_sdf'

# create a logger
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


def mt_to_vcf(
    batch: Batch,
    input_mt: str,
    samples: set[str],
    output_root: str,
) -> dict[str, tuple[str, hb.batch.job.Job | None]]:
    """
    takes a MT and converts to VCF
    optionally adds in extra header lines for VQSR filters
    returns a dictionary of all samples, their VCF paths
    if the VCF didn't already exist, this also contains a batch job

    Parameters
    ----------
    batch : batch to add jobs into
    input_mt : path to the MT to read into VCF
    samples : set of CPG sample IDs
    output_root :
    """

    # open the joint-call and check for the samples present
    init_batch()
    all_jc_samples = hl.read_matrix_table(input_mt).s.collect()
    samples_in_jc = set(all_jc_samples).intersection(samples)
    logging.info(f'Extracting {" ".join(samples_in_jc)} from the joint-call')
    sample_jobs = {}

    # extract all common samples into a separate file
    for sample in samples_in_jc:

        sample_path = os.path.join(
            output_root, 'single_sample_vcfs', f'{sample}.vcf.bgz'
        )

        if AnyPath(sample_path).exists():
            sample_jobs[sample] = (sample_path, None)
            logging.info(f'No action taken, {sample_path} already exists')
            continue

        job = batch.new_job(f'Extract {sample} from VCF')
        job.image(get_config()['workflow']['driver_image'])
        authenticate_cloud_credentials_in_job(job)
        prepare_git_job(
            job=job,
            organisation=get_organisation_name_from_current_directory(),
            repo_name=get_repo_name_from_current_directory(),
            commit=get_git_commit_ref_of_current_repository(),
        )

        job.command(
            f'PYTHONPATH=$(pwd) python3 {MT_TO_VCF_SCRIPT} '
            f'-i {input_mt} '
            f'-s {sample} '
            f'-o {sample_path} '
        )
        sample_jobs[sample] = (sample_path, job)

    return sample_jobs


def comparison_job(
    batch: Batch,
    dependency: hb.batch.job.Job | None,
    ss_vcf: str,
    sample: str,
    truth_vcf: str,
    truth_bed: str,
    comparison_folder: str,
    stratification: str | None = None,
):
    """

    Parameters
    ----------
    batch : the batch to run this job in
    dependency : a previous job to wait for
    ss_vcf : single sample VCF to validate
    sample : CPG ID
    truth_vcf : truth variant source
    truth_bed : confident region source
    comparison_folder :
    stratification : path to stratification BED data
    """

    job = batch.new_job(name=f'Compare {sample}')

    if dependency is not None:
        job.depends_on(dependency)

    job.image(image_path('happy'))
    job.memory('20Gi')
    job.storage('40Gi')
    job.cpu(2)
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

    # sdf loading as a Glob operation
    sdf = batch.read_input_group(
        **{file.name: file.as_uri() for file in AnyPath(REF_SDF).glob('*')}
    )

    # hap.py outputs:
    # output.extended.csv
    # output.vcf.gz
    # output.vcf.gz.tbi
    # output.roc.all.csv.gz
    # output.metrics.json.gz
    # output.runinfo.json
    # output.summary.csv

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
    command = (
        f'mkdir {job.output} && '
        f'hap.py {truth_input["vcf"]} {vcf_input["vcf"]} '
        f'-r {batch_ref["fasta"]} -R {truth_bed} '
        f'-o {job.output}/output --leftshift '
        f'--threads 4 --preprocess-truth '
        f'--engine-vcfeval-path=/opt/hap.py/libexec/rtg-tools-install/rtg '
        f'--engine-vcfeval-template {sdf} --engine=vcfeval '
    )

    # allow for stratification
    if stratification:
        strat_folder = to_path(stratification)
        assert strat_folder.exists(), (
            f'provided folder {stratification} does ' 'not exist, or was not accessible'
        )

        definitions = strat_folder / 'definition.tsv'
        assert definitions.exists(), (
            f'the region definition file ' f'{str(definitions)} does not exist'
        )

        strat_bed_files = list(strat_folder.glob('*.bed*'))
        assert (
            len(strat_bed_files) > 0
        ), 'There were no bed files in the stratified BED folder'

        # create a dictionary to pass to a the input generation
        strat_dict = {'definition.tsv': str(definitions)}
        strat_dict.update({file.name: str(file) for file in strat_bed_files})
        batch_beds = batch.read_input_group(**strat_dict)
        command += f'--stratification {batch_beds["definition.tsv"]}'

    job.command(command)

    # extract the results files
    batch.write_output(job.output, os.path.join(comparison_folder, sample))
    return job


def get_validation_samples() -> dict[str, str]:
    """
    query metamist for all sample IDs in this dataset

    Returns
    -------
    dict of {external: CPG ID}
    """

    party = ParticipantApi()
    results = party.get_external_participant_id_to_internal_sample_id(
        project=get_config()['workflow']['dataset']
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
        projects=[get_config()['workflow']['dataset']],
        sample_ids=[cpg_id],
        type=AnalysisType('custom'),
        meta={'type': 'validation'},
        active=True,
    )
    analyses = AnalysisApi().query_analyses(analysis_query_model=a_query_model)
    if len(analyses) > 1:
        raise Exception(
            f'Multiple [custom] analysis objects were found for '
            f'{cpg_id}, please set old analyses to active=False'
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


def post_results_job(
    batch: Batch,
    dependency: hb.batch.job.Job,
    sample_id: str,
    ss_vcf: str,
    truth_vcf: str,
    truth_bed: str,
    joint_mt: str,
    comparison_folder: str,
    stratified: str | None = None,
):
    """
    post the results to THE METAMIST using companion script

    Parameters
    ----------
    batch : batch to run the job in
    dependency : the analysis job we're dependent on
    sample_id : CPG ID
    ss_vcf : single sample VCF file used/created
    truth_vcf : source of truth variants
    truth_bed : source of confident regions
    joint_mt : joint-call MatrixTable
    comparison_folder :
    stratified : stratification files, if used
    """

    post_job = batch.new_job(name=f'Update metamist for {sample_id}')

    prepare_git_job(
        job=post_job,
        organisation=get_organisation_name_from_current_directory(),
        repo_name=get_repo_name_from_current_directory(),
        commit=get_git_commit_ref_of_current_repository(),
    )
    post_job.depends_on(dependency)
    post_job.image(get_config()['workflow']['driver_image'])
    copy_common_env(post_job)
    authenticate_cloud_credentials_in_job(post_job)
    job_cmd = (
        f'PYTHONPATH=$(pwd) python3 {RESULTS_SCRIPT} '
        f'--id {sample_id} '
        f'--folder {comparison_folder} '
        f'--ss {ss_vcf} '
        f'-t {truth_vcf} '
        f'-b {truth_bed} '
        f'--mt {joint_mt} '
    )

    # add stratification files if appropraite
    if stratified:
        job_cmd += f'--stratified {stratified}'

    post_job.command(job_cmd)


def main(input_file: str, stratification: str | None):
    """
    Parameters
    ----------
    input_file : path to the MT representing this joint-call
    stratification : the path to the stratification BED files
    """

    input_path = Path(input_file)
    if input_path.suffix != '.mt':
        logger.error('Expected a MT as input file')
        raise Exception('Expected a MT as input file')

    # # set the path for this output
    # process the MT to get the name
    validation_output_path = output_path(
        f'{input_path.name.replace(input_path.suffix, "")}'
    )

    validation_lookup = get_validation_samples()

    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = Batch(
        name='run validation bits and pieces',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_memory='highmem',
    )

    # generate single sample vcfs
    sample_jobs = mt_to_vcf(
        batch=batch,
        input_mt=input_file,
        samples=set(validation_lookup.keys()),
        output_root=validation_output_path,
    )

    if len(sample_jobs) == 0:
        raise Exception('No jobs/VCFs were created from this joint call')

    comparison_folder = os.path.join(validation_output_path, 'comparison')

    # iterate over the samples, and corresponding file paths/batch jobs
    for cpg_id, sample_data in sample_jobs.items():
        sample_vcf, vcf_job = sample_data

        # for each sample, use metamist to pull the corresponding truth and VCF
        # if there aren't any scheduled jobs, don't run the batch
        # more elegant solution? Like... a batch without jobs shouldn't try and run
        # or... batch.job_count()?
        truth_vcf, truth_bed = get_sample_truth(cpg_id)

        # we already assert that both are populated, so check one
        if truth_vcf is None:
            logger.error(f'Truth missing, skipping validation run for {cpg_id}')
            continue

        comparison = comparison_job(
            batch=batch,
            dependency=vcf_job,
            ss_vcf=sample_vcf,
            sample=cpg_id,
            truth_vcf=truth_vcf,
            truth_bed=truth_bed,
            comparison_folder=comparison_folder,
            stratification=stratification,
        )
        post_results_job(
            batch=batch,
            dependency=comparison,
            sample_id=cpg_id,
            ss_vcf=sample_vcf,
            truth_vcf=truth_vcf,
            truth_bed=truth_bed,
            joint_mt=input_file,
            comparison_folder=comparison_folder,
            stratified=stratification,
        )

    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('-i', help='input_path')
    parser.add_argument('-s', help='stratification BED directory')
    args = parser.parse_args()
    main(input_file=args.i, stratification=args.s)
