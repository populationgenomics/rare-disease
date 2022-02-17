#!/usr/bin/env python3


"""
wraps a slivar command submission
"""


import os
import hailtop.batch as hb

import click


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
SLIVAR_TAG = 'slivar:v0.2.7'
FULL_IMAGE = f'{AR_REPO}/{SLIVAR_TAG}'

# used to provide all VEP105 consequences, silences Slivar errors
CUSTOM_SLIVAR_CONS = 'gs://cpg-acute-care-test/vep/custom_slivar_order.txt'


@click.command()
@click.option(
    '--vcf',
    'vcf',
    help='vcf file to run the command on'
)
def main(vcf: str):
    """
    Create a Hail Batch, and run a tabix task within a job
    """

    vcf_out = os.path.splitext(vcf)[0] + '-slivar.vcf.bgz'

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(
        name='run_slivar',
        backend=service_backend
    )

    job = batch.new_job(name='run slivar')
    input_vcf_resource = batch.read_input(vcf)
    csq_list_resource = batch.read_input(CUSTOM_SLIVAR_CONS)

    job.command(
        f"""
        SLIVAR_IMPACTFUL_ORDER={csq_list_resource} slivar expr \
        --vcf {input_vcf_resource} \
        --pass-only \
        --skip-non-variable \
        --info 'INFO.impactful && variant.ALT[0] != "*"' \
        | bgzip -c -@ 4 > {job.out_vcf}
        """
    )
    job.cpu(2)
    job.memory('standard')  # ~ 4G/core ~ 7.5G
    job.storage('20G')
    batch.write_output(job.out_vcf, vcf_out)
    job.image(FULL_IMAGE)

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
