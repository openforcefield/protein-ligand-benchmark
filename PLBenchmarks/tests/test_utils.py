"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import utils
import pytest
from simtk import unit
import urllib
import sys

def test_findPdbUrl():
    """tests the findPdbUrl function"""
    pdbs = ['6R8L', '6N5O', '6SVL', '6J5R', '6J6D']
    for pdb in pdbs:
        assert f'REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3' == utils.findPdbUrl(pdb)

    strings = [ f'REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3' for pdb in pdbs ]
    for string in strings:
        assert string in utils.findPdbUrl(' '.join(pdbs)).split()

    with pytest.raises(ValueError, match=f'PDB fakepdb not found'):
        utils.findPdbUrl('fakepdb')

    pdbs[3] = 'fakepdb2'
    with pytest.raises(ValueError):
        utils.findPdbUrl(' '.join(pdbs))

def test_findDoiUrl():
    assert 'REP1http://dx.doi.org/10.1021/acs.jctc.8b00640REP2Mobley et al., J. Chem. Theory Comput. 2018REP3' == utils.findDoiUrl('10.1021/acs.jctc.8b00640')

    assert 'fakeDOI' == utils.findDoiUrl('fakeDOI')


def test_convertValue():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'dg'
    ##############################################    
    dg = unit.Quantity(1, unit.kilojoules_per_mole)
    assert pytest.approx(1.0, eps) == utils.convertValue(dg, 'dg', 'dg', outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(1.0, eps)  == utils.convertValue(dg, 'dg', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(0.239, eps) == utils.convertValue(dg, 'dg', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(0.6697, eps) == utils.convertValue(dg, 'dg', 'ki').value_in_unit(unit.molar)
    assert pytest.approx(0.6697, eps) == utils.convertValue(dg, 'dg', 'ic50').value_in_unit(unit.molar)
    assert pytest.approx(0.1741, eps) == utils.convertValue(dg, 'dg', 'pic50').value_in_unit(unit.dimensionless)
    with pytest.raises(NotImplementedError):
        assert '0.24 kcal/mol' == utils.convertValue(dg, 'dg', 'fakeObs')

    ##############################################
    # ORIGINAL = 'ki'
    ##############################################    
    ki = unit.Quantity(1, unit.molar)
    assert pytest.approx(0.0, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=300, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(0.0, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(0.0, eps) == utils.convertValue(ki, 'ki', 'dg', outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(1.0, eps) == utils.convertValue(ki, 'ki', 'ki', outUnit=unit.molar).value_in_unit(unit.molar)
    assert pytest.approx(1.0, eps) == utils.convertValue(ki, 'ki', 'ic50').value_in_unit(unit.molar)

    assert pytest.approx(0.0, eps) == utils.convertValue(ki, 'ki', 'pic50').value_in_unit(unit.dimensionless)

    with pytest.raises(NotImplementedError):
        assert 'xxx' == utils.convertValue(ki, 'ki', 'fakeObs')

    ki = unit.Quantity(1, unit.nano * unit.molar)
    assert pytest.approx(-51.69, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=300, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)

    assert pytest.approx(-12.35, eps) == utils.convertValue(ki, 'ki', 'dg').value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(-12.35, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=300).value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(-12.35, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=300, outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilocalories_per_mole)

    assert pytest.approx(-47.04, eps) == utils.convertValue(ki, 'ki', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)

    assert pytest.approx(1.00, eps) == utils.convertValue(ki, 'ki', 'ki').value_in_unit(unit.nano*unit.molar)
    assert pytest.approx(1.00, eps) == utils.convertValue(ki, 'ki', 'ki', outUnit=unit.nano*unit.molar).value_in_unit(unit.nano*unit.molar)

    assert pytest.approx(1.0, eps) == utils.convertValue(ki, 'ki', 'ic50').value_in_unit(unit.nano * unit.molar)

    assert pytest.approx(9, eps) == utils.convertValue(ki, 'ki', 'pic50').value_in_unit(unit.dimensionless)

    with pytest.raises(NotImplementedError):
        assert 'xxx' == utils.convertValue(ki, 'ki', 'fakeObs')



    ##############################################
    # ORIGINAL = 'ic50'
    ##############################################    
    ic50 = unit.Quantity(1, unit.molar)
    assert pytest.approx(0.0, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=300, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(0.0, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(0.0, eps) == utils.convertValue(ic50, 'ic50', 'dg', outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilojoules_per_mole)

    assert pytest.approx(1.0, eps) == utils.convertValue(ic50, 'ic50', 'ki', outUnit=unit.molar).value_in_unit(unit.molar)

    assert pytest.approx(1.0, eps) == utils.convertValue(ic50, 'ic50', 'ic50', outUnit=unit.nano*unit.molar).value_in_unit(unit.molar)

    assert pytest.approx(0.0, eps) == utils.convertValue(ic50, 'ic50', 'pic50').value_in_unit(unit.dimensionless)

    with pytest.raises(NotImplementedError):
        utils.convertValue(ic50, 'ic50', 'fakeObs')


    ic50 = unit.Quantity(1, unit.nano * unit.molar)
    assert pytest.approx(-51.69, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=300, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)
    assert pytest.approx(-12.35, eps) == utils.convertValue(ic50, 'ic50', 'dg').value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(-12.35, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=300).value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(-12.35, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=300, outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilocalories_per_mole)
    assert pytest.approx(-47.04, eps) == utils.convertValue(ic50, 'ic50', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole)

    assert pytest.approx(1.00, eps) == utils.convertValue(ic50, 'ic50', 'ki').value_in_unit(unit.nano*unit.molar)
    assert pytest.approx(1.00, eps) == utils.convertValue(ic50, 'ic50', 'ki', outUnit=unit.nano*unit.molar).value_in_unit(unit.nano*unit.molar)

    assert pytest.approx(1.00, eps) == utils.convertValue(ic50, 'ic50', 'ic50').value_in_unit(unit.nano*unit.molar)
    assert pytest.approx(1.00, eps) == utils.convertValue(ic50, 'ic50', 'ic50', outUnit=unit.nano*unit.molar).value_in_unit(unit.nano*unit.molar)

    assert pytest.approx(9.00, eps) == utils.convertValue(ic50, 'ic50', 'pic50').value_in_unit(unit.dimensionless)

    with pytest.raises(NotImplementedError):
        utils.convertValue(ic50, 'ic50', 'fakeObs')



    ##############################################
    # ORIGINAL = 'pic50'
    ##############################################    
    pic50 = unit.Quantity(0, unit.dimensionless)
    assert pytest.approx(0.0, eps)   == utils.convertValue(pic50, 'pic50', 'dg', temperature = 300.0, outUnit = unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole) 
    assert pytest.approx(0.0, eps)   == utils.convertValue(pic50, 'pic50', 'dg', temperature = 273.0, outUnit = unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole) 
    assert pytest.approx(0.0, eps)   == utils.convertValue(pic50, 'pic50', 'dg', outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilojoules_per_mole) 

    assert pytest.approx(1.0, eps)   == utils.convertValue(pic50, 'pic50', 'ki', outUnit=unit.molar).value_in_unit(unit.molar) 

    assert pytest.approx(1.0, eps)   == utils.convertValue(pic50, 'pic50', 'ic50', outUnit=unit.molar).value_in_unit(unit.molar) 

    assert pytest.approx(0.0, eps)      == utils.convertValue(pic50, 'pic50', 'pic50').value_in_unit(unit.dimensionless)

    with pytest.raises(NotImplementedError):
        utils.convertValue(pic50, 'pic50', 'fakeObs')

    pic50 = unit.Quantity(9, unit.dimensionless)
    assert pytest.approx(-51.69, eps) == utils.convertValue(pic50, 'pic50', 'dg', temperature=300, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole) 

    assert pytest.approx(-12.35, eps) == utils.convertValue(pic50, 'pic50', 'dg').value_in_unit(unit.kilocalories_per_mole) 
    assert pytest.approx(-12.35, eps) == utils.convertValue(pic50, 'pic50', 'dg', temperature=300).value_in_unit(unit.kilocalories_per_mole) 
    assert pytest.approx(-12.35, eps) == utils.convertValue(pic50, 'pic50', 'dg', temperature=300, outUnit=unit.kilocalories_per_mole).value_in_unit(unit.kilocalories_per_mole) 

    assert pytest.approx(-47.04, eps) == utils.convertValue(pic50, 'pic50', 'dg', temperature=273, outUnit=unit.kilojoules_per_mole).value_in_unit(unit.kilojoules_per_mole) 

    assert pytest.approx(1.00, eps)   == utils.convertValue(pic50, 'pic50', 'ki').value_in_unit(unit.nano*unit.molar) 

    assert pytest.approx(1.00, eps)   == utils.convertValue(pic50, 'pic50', 'ic50').value_in_unit(unit.nano * unit.molar) 
    assert pytest.approx(1.00, eps)   == utils.convertValue(pic50, 'pic50', 'ic50', outUnit=unit.nano*unit.molar).value_in_unit(unit.nano * unit.molar) 

    assert pytest.approx(9, eps)      == utils.convertValue(pic50, 'pic50', 'pic50').value_in_unit(unit.dimensionless) 
    with pytest.raises(NotImplementedError):
        utils.convertValue(pic50, 'pic50', 'fakeObs')
