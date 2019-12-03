"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
import PLBenchmarks
import pytest
import sys

def test_PLBenchmarks_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "PLBenchmarks" in sys.modules
