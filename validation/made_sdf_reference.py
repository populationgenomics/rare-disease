#!/usr/bin/env python3

"""
wraps the validation script(s)
"""

import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, remote_tmpdir, image_path

REF_GENOME = (
    'gs://cpg-reference/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta'
)


if __name__ == '__main__':

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
    job = batch.new_job(name='Generate Reference SDF')
    job.image(image_path('happy-vcfeval'))
    job.memory('20Gi')
    batch_ref = batch.read_input_group(
        **{'fasta': REF_GENOME, 'index': f'{REF_GENOME}.fai'}
    )

    job.declare_resource_group(
        output={
            'format.log': '{root}/format.log',
            'mainIndex': '{root}/mainIndex',
            'namedata0': '{root}/namedata0',
            'nameIndex0': '{root}/nameIndex0',
            'namepointer0': '{root}/namepointer0',
            'reference.txt': '{root}/reference.txt',
            'seqdata0': '{root}/seqdata0',
            'seqdata1': '{root}/seqdata1',
            'seqdata2': '{root}/seqdata2',
            'seqdata3': '{root}/seqdata3',
            'seqpointer0': '{root}/seqpointer0',
            'seqpointer1': '{root}/seqpointer1',
            'seqpointer2': '{root}/seqpointer2',
            'seqpointer3': '{root}/seqpointer3',
            'sequenceIndex0': '{root}/sequenceIndex0',
            'progress': '{root}/progress',
            'summary.txt': '{root}/summary.txt',
        }
    )

    # in future don't regenerate SDF...
    job.command(
        'java -jar -Xmx16G /vcfeval/RTG.jar '
        f'format -o {job.output} {batch_ref["fasta"]}'
    )
    batch.write_output(job.output, output_path('ref_genome_SDF'))
    batch.run(wait=False)
