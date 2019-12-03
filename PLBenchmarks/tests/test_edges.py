"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import edges
from simtk import unit
import pytest
import pandas as pd
import yaml
try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text


def test_edges():
    df = edges.getEdgesSet('jnk1')
    assert df.shape[0] > 0
    assert df.shape[1] == 2
