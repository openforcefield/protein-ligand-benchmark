"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import util
import pytest
from simtk import unit
import urllib
import sys

def test_findPdbUrl():
    """tests the findPdbUrl function"""
    pdbs = ['6R8L', '6N5O', '6SVL', '6J5R', '6J6D']
    for pdb in pdbs:
        assert f'REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3' == util.findPdbUrl(pdb)

    strings = [ f'REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3' for pdb in pdbs ]
    for string in strings:
        assert string in util.findPdbUrl(' '.join(pdbs)).split()

    with pytest.raises(ValueError, match=f'PDB fakepdb not found'):
        util.findPdbUrl('fakepdb')

    pdbs[3] = 'fakepdb2'
    with pytest.raises(ValueError):
        util.findPdbUrl(' '.join(pdbs))

def test_findDoiUrl():
    assert 'REP1http://dx.doi.org/10.1021/acs.jctc.8b00640REP2Mobley et al., J. Chem. Theory Comput. 2018REP3' == util.findDoiUrl('10.1021/acs.jctc.8b00640')

    assert 'fakeDOI' == util.findDoiUrl('fakeDOI')


def test_convertValue():
    # original = 'dg'
    dg = unit.Quantity(1, unit.kilojoules_per_mole)
    assert '1.00 kJ/mol' == util.convertValue(dg, 'dg', 'dg', energyUnit=unit.kilojoules_per_mole)
    assert '0.24 kcal/mol' == util.convertValue(dg, 'dg', 'dg', energyUnit=unit.kilocalories_per_mole)

    with pytest.raises(NotImplementedError):
        assert '1.00 kJ/mol' == util.convertValue(dg, 'dg', 'ki')
    with pytest.raises(NotImplementedError):
        assert '0.24 kcal/mol' == util.convertValue(dg, 'dg', 'ic50')
    with pytest.raises(NotImplementedError):
        assert '0.24 kcal/mol' == util.convertValue(dg, 'dg', 'pic50')
    with pytest.raises(NotImplementedError):
        assert '0.24 kcal/mol' == util.convertValue(dg, 'dg', 'fakeObs')