"""
test class for the utils collection
"""

import os
import pytest

from cyvcf2 import VCFReader
from reanalysis.utils import AnalysisVariant

PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')
PANELAPP_LATEST = os.path.join(INPUT, 'panelapp_current_137.json')


@pytest.fixture(name='variant_list')
def fixture_variant_list():
    """

    :return:
    """
