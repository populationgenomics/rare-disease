"""
a pipe-thru program that can filter lines based on panelapp dates
requires static data, so it partners with parse_panelapp_data.py

bcftools view MY_VCF.vcf.bgz <region> |
    python date_check.py
    input/panelapp_137_dates.json
    2021-01-01
> newfile.vcf

removes rows containing variants in genes added to
the specified panel prior to the datetime given

genes without entries in the JSON will default to the threshold
i.e. will pass into the output

Header lines ignored

Needs to work in python2 inside the container
the container doesn't have click or argparse installed...

this only works for the Hail formatted VCF export
"""

import json
import re
import sys
from datetime import datetime


ENSG_RE = re.compile(r'vep_gene_id=(?P<gene>ENSG[0-9]*);')


def read_pipe():
    """
    generator for lines coming through STDIN
    :return:
    """
    try:
        for line in sys.stdin:
            if line.startswith("#"):
                print(line.rstrip())
                continue
            yield line
    except KeyboardInterrupt:
        sys.stdout.flush()


def main(
    panel_json='/data_in/input/panelapp_137_dates.json', date_threshold='2021-01-01'
):
    """

    :param panel_json:
    :param date_threshold:
    :return:
    """
    threshold = datetime.strptime(date_threshold, '%Y-%m-%d')

    # read input JSON, saving {gene: datetime} lookup
    dates = {}
    with open(panel_json, 'r') as handle:  # pylint: disable=W1514
        date_dict = json.load(handle)
        # wrote this a dictionary comprehension and python2 didn't like it
        for key in date_dict.keys():
            dates[key] = datetime.strptime(date_dict[key]['since'], '%Y-%m-%d')

    for line in read_pipe():
        gene = ENSG_RE.search(line).group('gene')
        if dates.get(gene, threshold) >= threshold:
            print(line.rstrip())


if __name__ == "__main__":
    main(panel_json=sys.argv[1], date_threshold=sys.argv[2])
