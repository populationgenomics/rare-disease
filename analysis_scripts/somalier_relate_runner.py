# pylint: disable=missing-function-docstring,no-member
"""
** Adapted from Hope's script **
https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/sex_ploidy/somalier_relate_runner.py

This script runs somalier relate on a set of cram.somalier files provided as input directory file path.
Alternately, you can specify individual somalier files to compare.
# TODO: Allow specifying different inputs (e.g. Sequencing Group IDs, external Family IDs).

Comparing two directories:

```bash
 analysis-runner --dataset "dataset1" \
    --description "Somalier relate runner" \
    --access-level "standard" \
    --output-dir "qc-stand-alone/dataset1/somalier" \
    somalier_relate_runner.py --input-dir-1=gs://cpg-dataset1-main/cram --input-dir-2=gs://cpg-dataset2-main/cram
```

Comparing two files:

```bash
 analysis-runner --dataset "dataset1" \
    --description "Somalier relate runner" \
    --access-level "standard" \
    --output-dir "qc-stand-alone/dataset1/somalier" \
    somalier_relate_runner.py -i gs://cpg-dataset1-main/cram/CRAM1.cram.somalier -i gs://cpg-dataset2-main/cram/CRAM2.cram.somalier
```
"""

import click
from cpg_utils import to_path
from cpg_utils.config import config_retrieve, output_path
from cpg_utils.hail_batch import get_batch

SOMALIER_IMAGE = config_retrieve(['images', 'somalier'])


@click.option(
    '--input-dir-1',
    help='Input directory to cram.somalier files',
    default=None,
)
@click.option(
    '--input-dir-2',
    help='2nd (optional)input directory to cram.somalier files',
    default=None,
)
@click.option('-p', '--expected-ped-path', help='Expected pedigree file', default=None)
@click.option(
    '-i',
    '--input-filepaths',
    help='Input path to cram.somalier files',
    multiple=True,
)
@click.option(
    '--job-storage',
    help='Storage of the Hail batch job eg 30G',
    default='10G',
)
@click.option('--job-memory', help='Memory of the Hail batch job', default='8G')
@click.option('--job-ncpu', help='Number of CPUs of the Hail batch job', default=4)
@click.command()
def main(
    job_memory: str,
    job_ncpu: int,
    job_storage: str,
    input_filepaths: str,
    expected_ped_path: str | None,
    input_dir_1: str | None,
    input_dir_2: str | None,
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    if not input_dir_1 and not input_dir_2:
        input_files = list(input_filepaths)
        if any(not file.endswith('.somalier') for file in input_files):
            raise ValueError('Please provide .somalier files only')
        if len(input_files) < 2:  # noqa: PLR2004
            raise ValueError('Please provide at least 2 input files for comparison')
    else:
        input_files = list(to_path(input_dir_1).glob('*.somalier'))
        input_files = [
            str(gs_path) for gs_path in input_files
        ]  # coverts into a string type
        if input_dir_2 is not None:
            input_files_2 = list(to_path(input_dir_2).glob('*.somalier'))
            input_files_2 = [str(gs_path) for gs_path in input_files_2]
            input_files.extend(input_files_2)

    num_samples = len(input_files)
    batch_input_files = []
    for each_file in input_files:
        batch_input_files.append(b.read_input(each_file))

    if expected_ped_path:
        expected_ped_path = b.read_input(expected_ped_path)

    somalier_job = b.new_job(name=f'Somalier relate: {num_samples} samples')
    somalier_job.image(SOMALIER_IMAGE)
    somalier_job.storage(job_storage)
    somalier_job.memory(job_memory)
    somalier_job.cpu(job_ncpu)

    somalier_job.declare_resource_group(
        output={
            'pairs': '{root}.pairs.tsv',
            'samples': '{root}.samples.tsv',
            'html': '{root}.html',
        },
    )

    expected_ped_path = f'--ped {expected_ped_path}' if expected_ped_path else ''
    somalier_job.command(
        f"""
            somalier relate  \\
            {" ".join(batch_input_files)} \\
            {expected_ped_path} --infer \\
            -o {somalier_job.output}
        """,
    )

    # Output writing
    out_path = output_path(f'{num_samples}_samples_somalier', 'analysis')
    b.write_output(somalier_job.output, out_path)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
