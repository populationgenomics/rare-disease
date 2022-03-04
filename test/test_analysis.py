"""
test class for the class validation class
"""
import os
import pytest
import cyvcf2

from reanalysis.validate_classifications import (
    parse_comp_hets,
    read_json_dictionary,
)
from reanalysis.utils import (
    string_format_variant,
)

PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
JSON_STUB = os.path.join(INPUT, 'json_example.json')
COMP_HET = os.path.join(INPUT, 'test_comp_het.vcf.bgz')
SINGLE_VAR = os.path.join(INPUT, 'single_var.vcf.bgz')

# variants within the compound het file
# string value will be equivalent to the repr() of cyvcf2.Variant
expected_output = {
    'CPG54445': {
        '8-143919896-G-A': 'Variant(chr8:143925556 C/T)',
        '8-143925556-C-T': 'Variant(chr8:143919896 G/A)',
    },
    'CPG54569': {
        '2-177663913-G-A': 'Variant(chr2:178071545 T/C)',
        '2-178071545-T-C': 'Variant(chr2:177663913 G/A)',
    },
    'CPG11510': {
        '7-139875609-C-T': 'Variant(chr7:140020029 C/T)',
        '7-140020029-C-T': 'Variant(chr7:139875609 C/T)',
    },
}


def test_read_json():
    """

    :return:
    """
    assert read_json_dictionary(JSON_STUB) == {'key': 'value'}
    with pytest.raises(AssertionError):
        read_json_dictionary('not_a_real_path')


def test_comp_het_gather():
    """
    check that we gather the right data from the comphet het vcf
    :return:
    """
    result = parse_comp_hets(COMP_HET)
    assert set(expected_output.keys()) == set(result.keys())
    for sample, variants in result.items():
        assert set(expected_output[sample].keys()) == set(variants.keys())
        for var_key, variant in variants.items():
            assert variant.string == expected_output[sample][var_key]


def test_variant_string_format():
    """
    open a single variant test VCF
    chr2	177663913	rs61306957	G	A
    :return:
    """
    for variant in cyvcf2.VCFReader(SINGLE_VAR):
        assert string_format_variant(variant) == '2-177663913-G-A'
        assert (
            string_format_variant(variant, transcript=True)
            == '2-177663913-G-A-ENST00000286063'
        )
        break
