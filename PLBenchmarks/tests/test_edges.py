"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import edges



def test_edgeSet():
    edg = edges.edgeSet("mcl1")
    df = edg.getDF()
    html = edg.getHTML()
    d = edg.getDict()
