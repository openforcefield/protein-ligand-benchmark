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


def test_edgeSet():
    edg = edges.edgeSet('jnk1')
    df = edg.getDF()
    html = edg.getHTML()
    d = edg.getDict()

