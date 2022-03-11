"""
test class for the class validation class
"""
import json
import os

from reanalysis.validate_classifications import parse_comp_hets


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
JSON_STUB = os.path.join(INPUT, 'json_example.json')
COMP_HET = os.path.join(INPUT, 'comp_het.vcf.bgz')
SINGLE_VAR = os.path.join(INPUT, 'single_hail.vcf.bgz')

with open(os.path.join(INPUT, 'config.json'), 'r') as handle:
    FULL_CONF = json.load(handle)

# variants within the compound het file
# string value will be equivalent to the repr() of cyvcf2.Variant
expected_output = {
    'CPG54445': {
        '8-143919896-G-A': '8-143925556-C-T',
        '8-143925556-C-T': '8-143919896-G-A',
    },
    'CPG54569': {
        '2-177663913-G-A': '2-178071545-T-C',
        '2-178071545-T-C': '2-177663913-G-A',
    },
    'CPG11791': {
        '7-143342478-A-G': '7-143330817-G-A',
        '7-143330817-G-A': '7-143342478-A-G',
    },
}


def test_comp_het_gather():
    """
    check that we gather the right data from the comphet het vcf
    :return:
    """
    result = parse_comp_hets(COMP_HET, FULL_CONF)
    assert set(expected_output.keys()) == set(result.keys())
    for sample, variants in result.items():
        assert set(expected_output[sample].keys()) == set(variants.keys())
        for var_key, variant in variants.items():
            assert variant.string == expected_output[sample][var_key]
