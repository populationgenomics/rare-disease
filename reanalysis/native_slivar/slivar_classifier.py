#!/usr/bin/env python3


"""
wraps multiple slivar command submissions
one command launched per class
"""


from typing import Optional
import os

import click
import hailtop.batch as hb


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
SLIVAR_TAG = 'slivar:v0.2.7'
SLIVAR_IMAGE = f'{AR_REPO}/{SLIVAR_TAG}'
SLIVAR_TEMPLATE = """
set -ex
SLIVAR_QUIET="true" slivar expr \
--vcf {input_vcf} \
--pass-only \
--skip-non-variable \
--info '{info_expr}' \
--sample-expr 'het:sample.DP && sample.GQ > 10 && sample.het' \
--sample-expr 'hom:sample.DP && sample.GQ > 10 && sample.hom_alt' \
| bgzip -c -@ 4 > {output}
tabix {output}
"""

CLINVAR_PATHO = ' || '.join(
    [f'INFO.clinvar_sig == "{csq}' for csq in ["Pathogenic", "Likely_pathogenic"]]
)


# used to provide all VEP105 consequences, silences Slivar errors
# not required for Hail-prepared slivar
CUSTOM_SLIVAR_CONS = 'gs://cpg-acute-care-test/vep/custom_slivar_order.txt'


def make_generic_job(batch: hb.Batch, name: str):
    """
    returns a new, appropriately resourced job
    :param batch:
    :param name:
    :return:
    """
    job = batch.new_job(name=name)
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')
    job.image(SLIVAR_IMAGE)
    return job


def make_class_1_job(
    batch: hb.Batch, input_vcf: str, name: Optional[str] = 'run slivar class 1'
):
    """
    relatively simple,
    - Pathogenic or Likely_pathogenic
    - at least one clinvar star

    :param batch:
    :param input_vcf:
    :param name:
    :return:
    """
    job = make_generic_job(batch, name=name)

    job.declare_resource_group(
        vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    job_cmd = SLIVAR_TEMPLATE.format(
        input_vcf=input_vcf,
        info_expr=f'({CLINVAR_PATHO}) && INFO.clinvar_stars > 0',
        output=job.vcf["vcf.bgz"],
    )
    job.command(job_cmd)

    return job


def make_class_3_job(
    batch: hb.Batch, input_vcf: str, name: Optional[str] = 'run slivar class 3'
):
    """
    High impact consequence, rare, and supported in Clinvar (any quality)
    '(
        variant.INFO.vep_csq == "transcript_ablation" ||
        variant.INFO.vep_csq == "splice_acceptor" ||
        variant.INFO.vep_csq == "splice_donor" ||
        variant.INFO.vep_csq == "stop_gained" ||
        variant.INFO.vep_csq == "start_lost" ||
        variant.INFO.vep_csq == "frameshift" ||
        variant.INFO.vep_csq == "stop_lost"
    ) &&
    (INFO.gnomad_af < 0.002) &&
    (
        variant.INFO.clinvar_sig == "Pathogenic" ||
        variant.INFO.clinvar_sig == "Likely_pathogenic"
    )'

    :param batch:
    :param input_vcf:
    :param name:
    :return:
    """
    job = make_generic_job(batch, name=name)

    job.declare_resource_group(
        vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    all_csq_conditions = ' || '.join(
        [
            f'variant.INFO.vep_csq == "{csq}"'
            for csq in [
                'transcript_ablation',
                "splice_acceptor",
                "splice_donor",
                "stop_gained",
                "start_lost",
                "frameshift",
                "stop_lost",
            ]
        ]
    )

    job_cmd = SLIVAR_TEMPLATE.format(
        input_vcf=input_vcf,
        info_expr=f'({all_csq_conditions}) && (INFO.gnomad_af < 0.002) && ({CLINVAR_PATHO})',
        output=job.vcf["vcf.bgz"],
    )
    job.command(job_cmd)

    return job


def make_class_4_job(
    batch: hb.Batch, input_vcf: str, name: Optional[str] = 'run slivar class 4'
):
    """
    mostly about in silico prediction - either cadd, revel, or splice AI strong
    OR support across all the other tools (probably need more tools for confidence)
      (
      INFO.vep_sift_prediction == "deleterious" &&
      INFO.vep_sift_score < 0.04 &&
      INFO.vep_polyphen_score >= 0.89 &&
      INFO.vep_polyphen_prediction == "probably_damaging"
    ) || (
      INFO.cadd > 30 || INFO.revel > 0.4 || INFO.splice_ai_delta >= 0.5
    ) && (
      INFO.gnomad_af < 0.002
    )
      :param batch:
      :param input_vcf:
      :param name:
      :return:
    """
    job = make_generic_job(batch, name=name)

    job.declare_resource_group(
        vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )

    # requires consensus
    soft_toolgroup = (
        '('
        'INFO.vep_sift_prediction == "deleterious" && '
        'INFO.vep_sift_score < 0.04 && '
        'INFO.vep_polyphen_score >= 0.89 && '
        'INFO.vep_polyphen_prediction == "probably_damaging"'
        ')'
    )

    # requires match for one of these, we like these more
    hard_toolgroup = (
        '(INFO.cadd > 30 || INFO.revel > 0.4 || INFO.splice_ai_delta >= 0.5)'
    )

    job_cmd = SLIVAR_TEMPLATE.format(
        input_vcf=input_vcf,
        info_expr=f'INFO.gnomad_af < 0.002 && ({hard_toolgroup} || {soft_toolgroup})',
        output=job.vcf["vcf.bgz"],
    )
    job.command(job_cmd)

    return job


@click.command()
@click.option('--vcf', 'vcf', help='vcf file to run the command on')
@click.option('--bucket', 'bucket', help='gcp bucket to write output into')
def main(vcf: str, bucket: str):
    """
    Create a Hail Batch, and run slivar tasks in a batch

    This process assumes that the VCF file has been reduced in Hail
    e.g. instead of a CSQ field to parse, each variant should be
    reduced to a single consequence, under strict criteria
    and the individual csq elements should be exposed as info fields
    :param vcf: input VCF - Hail-VEP export expected
    :param bucket: GCP path to write outputs to
    :return:
    """

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        bucket=os.getenv('HAIL_BUCKET'),
    )

    # create a hail batch
    batch = hb.Batch(name='run_slivar', backend=service_backend)

    # read VCF in once for this batch
    vcf_resource = batch.read_input(vcf)

    # create the job & write out the full resource group, including the index
    for method, name in [
        (make_class_1_job, 'slivar_class_1'),
        (make_class_3_job, 'slivar_class_3'),
        (make_class_4_job, 'slivar_class_4'),
    ]:
        class_job = method(batch=batch, input_vcf=vcf_resource, name=name)
        batch.write_output(class_job.vcf, os.path.join(bucket, name))

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
