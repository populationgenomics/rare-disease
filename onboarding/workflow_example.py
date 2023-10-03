#!/usr/bin/env python3


"""
A skeleton of a Hail Batch workflow, with examples of a few core concepts

General design -
1. Make a Metamist query to find a group of samples
2. Make a metamist query to find some files assc. with those samples
    n.b. this could be one step, but it's easier to document
3. Start a Hail Batch
    - Specify a docker image to use
    - Specify any other parameters (memory, cpu, etc)
4. Add jobs to the batch, executing some command on files/groups of files
5. Write the output of those jobs to GCP

This is on the basis that we are running this using AnalysisRunner, so
aspects of the job set up (e.g. setting a project, providing in/out paths)
are at least partially handled by a range of helper methods
"""

from argparse import ArgumentParser

# from metamist.graphql import gql, query
# from collections import defaultdict
# from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, output_path

from cpg_workflows.batch import get_batch


def get_samples(samples: list[str]) -> list[str]:
    """
    go to metamist, find some samples
    Args:
        samples ():

    Returns:

    """

    # project = get_config()['workflow']['dataset']
    return samples


def get_files(samples: list[str]) -> dict[str, list[str]]:
    """
    go to metamist, find some files
    Args:
        samples (list[str]): CPG IDs

    Returns:
        Dict, keyed by sample ID, of lists of file paths
        {
            CPGAAAA: [
                gs://cpg-fewgenomes-test/analysis/2021-06-01/CPGAAAA_1.fq.gz,
                gs://cpg-fewgenomes-test/analysis/2021-06-01/CPGAAAA_2.fq.gz,
            ],
        }
    """

    # project = get_config()['workflow']['dataset']

    # do a metatmist query, or don't, I don't care
    ret_dict: dict[str, list[str]] = {}
    for sample in samples:
        ret_dict[sample] = []  # integrate query results

    return ret_dict


def main(samples: list[str]):
    """
    Control method for this process

    Returns:

    """
    # region: metadata queries
    # This section is all executed prior to the workflow being scheduled,
    # so we have access to all variables
    # Metamist Query for samples, or just take some from CLI args
    samples = get_samples(samples)

    # Metamist query for files
    file_dict = get_files(samples)
    # endregion

    # region: Batch time
    b = get_batch('Schedule some worker tasks')
    for sample, files in file_dict.items():
        # check that we have 2 files (assuming FQ, rather than BAM)
        assert len(files) == 2

        # Create a job for each sample
        j = b.new_job(f'Job for {sample}')

        # Set the docker image to use in this job
        # this pulls the image path from the portion of the config
        # populated by the images repository
        j.image(image_path('fastqe'))

        # set some other job attributes, if required
        # by default, ~4GB RAM, 0 additional storage, 1 CPU
        # j.cpu(2)
        # j.memory('4Gi')
        # j.storage('10Gi')

        # read data into the batch tmp resource location
        file_1 = b.read_input(files[0])
        file_2 = b.read_input(files[1])

        # Set the command to run
        # batch.read_input will create a new path like /io/batch/75264c/CPGAAAA_1.fq.gz
        # accessible from inside the job container, and unique to this batch/job
        j.command(
            f'echo "Hello world, I am a job for {sample}!, using {file_1} and {file_2}"'
            f'I\'m also creating an output file at {j.output}'
            f'echo "Some outputs" > {j.output}'
        )

        # read the output out into GCP
        # The helper method output_path() will create a new path based on the current project,
        # test/main, and the output prefix provided to analysis_runner
        # -o my_output
        # --dataset my-dataset
        # --access_level test
        # output_path('this_file.txt')
        # -> gs://cpg-my-dataset-test/my_output/this_file.txt
        b.write_output(j.output, output_path(f'/{sample}.txt'))
    # endregion


if __name__ == '__main__':
    # optional direct input of samples
    parser = ArgumentParser(description='Hail Batch workflow Skeleton')
    parser.add_argument('-s', nargs='+', help='Sample IDs to run on', required=True)
    args, fails = parser.parse_known_args()

    if fails:
        raise ValueError(f'Failed to parse arguments: {fails}')
    main(samples=args.s)
