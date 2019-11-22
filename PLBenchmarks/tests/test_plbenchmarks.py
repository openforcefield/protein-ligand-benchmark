"""
Unit and regression test for the plbenchmarks package.
"""

# Import package, test suite, and other packages as needed
import plbenchmarks
import pytest
import sys

def test_plbenchmarks_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "plbenchmarks" in sys.modules
