# !/usr/bin/env python3


"""
script to run slivar commands within a cloud instance
uses the default slivar container, so split-vep not present
"""


from typing import Optional
import os.path

import click
from hailtop.batch.job import Job

from cpg_pipes import hailbatch

SLIVAR_IMAGE = 'australia-southeast1-docker_slivar.pkg.dev/cpg-common/images/slivar:v0.2.7'

# file containing standard slivar csq ordering,
# plus all sub-consequences added in VEP 105
# this prevents errors being written for almost every variant
CUSTOM_SLIVAR_CONS = 'gs://cpg-acute-care-test/vep/custom_slivar_order.txt'

@click.command()
@click.option(
    '--vcf',
    'vcf_in',
    required=True,
)
@click.option(
    '--project',
    'project',
)
def main(vcf: str, project: str):
    """
    takes a VCF, and runs the enclosed Slivar command on it
    :param vcf:
    :param project:
    :return:
    """

    batch = hailbatch.setup_batch('Run Slivar', analysis_project_name=project)
    vcf_out = os.path.splitext(vcf)[0] + '-slivar.vcf.bgz'
    _job = slivar(batch=batch, vcf_path=vcf, out_vcf_path=vcf_out)
    batch.run(wait=False)


def slivar(
        batch,
        vcf_path: str,
        out_vcf_path: Optional[str] = None,
) -> Job:
    """
    Runs Slivar on provided VCF. Create a new VCF by altering input of extension
    Use a GSUtil copy in for input data
    """
    j = batch.new_job(f'Slivar_on_{vcf_path}')
    j.image(SLIVAR_IMAGE)
    hailbatch.STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    cmd = f"""\
    retry_gs_cp {vcf_path} input.vcf.gz
    retry_gs_cp {vcf_path}.tbi input.vcf.gz.tbi
    retry_gs_cp {CUSTOM_SLIVAR_CONS} custom_cons.txt
    SLIVAR_IMPACTFUL_ORDER=custom_cons.txt slivar expr \
    --vcf input.vcf.gz \
    --pass-only \
    --skip-non-variable \
    --info '
    (INFO.CLIN_SIG == "pathogenic" || INFO.CLIN_SIG == "likely_pathogenic")
    && variant.ALT[0] != "*"
    && (INFO.MAX_AF == "." || INFO.MAX_AF < 0.01)
    ' \
    | bgzip -c -@ 4 > {j.out_vcf}
    """
    j.command(
        hailbatch.wrap_command(
            cmd,
            monitor_space=True,
            define_retry_function=True
        )
    )
    if out_vcf_path:
        batch.write_output(j.out_vcf, out_vcf_path)
    return j


if __name__ == '__main__':
    # entrypoint
    main()  # pylint: disable=E1120
