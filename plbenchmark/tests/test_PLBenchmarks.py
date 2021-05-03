"""
Unit and regression test for the plbenchmark package.
"""

# Import package, test suite, and other packages as needed
import plbenchmark
import pytest
import sys


def test_plbenchmark_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "plbenchmark" in sys.modules
