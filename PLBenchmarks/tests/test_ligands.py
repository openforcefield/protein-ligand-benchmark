"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import ligands
from simtk import unit
import pytest
import pandas as pd
import yaml
try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text




def test_ligand():
    file = open_text('PLBenchmarks.data.01_jnk1.00_data', 'ligands.yml')    
    data = yaml.full_load_all(file)
    dfs=[]
    for d in data:
        l = ligands.ligand(d)
        l.deriveObservables(derivedObs='dg')
        l.addMolToFrame()
        dfs.append(l.getDF(['name', 'ROMol', 'DerivedMeasurement']))
    df = pd.DataFrame(dfs)
    assert df.shape[0] == 21

    for n in df.name:
        assert n in ['lig_17124-1',
                     'lig_18624-1',
                     'lig_18625-1',
                     'lig_18626-1',
                     'lig_18627-1',
                     'lig_18628-1',
                     'lig_18629-1',
                     'lig_18630-1',
                     'lig_18631-1',
                     'lig_18632-1',
                     'lig_18633-1',
                     'lig_18634-1',
                     'lig_18635-1',
                     'lig_18636-1',
                     'lig_18637-1',
                     'lig_18638-1',
                     'lig_18639-1',
                     'lig_18652-1',
                     'lig_18658-1',
                     'lig_18659-1',
                     'lig_18660-1']
        


    # Check whether the values in the repo are the same and correctly converted by comparing to the values in the JACS paper
    jacs_data = {'lig_18628-1' : -8.7,
                 'lig_18624-1' : -8.49,
                 'lig_18639-1' : -9.74,
                 'lig_18660-1' : -8.7,
                 'lig_18630-1' : -9.14,
                 'lig_18632-1' : -9.08,
                 'lig_18636-1' : -7.51,
                 'lig_18652-1' : -10.68,
                 'lig_17124-1' : -9.68,
                 'lig_18635-1' : -7.29,
                 'lig_18627-1' : -8.48,
                 'lig_18637-1' : -10.14,
                 'lig_18634-1' : -9.99,
                 'lig_18629-1' : -8.67,
                 'lig_18631-1' : -9.41,
                 'lig_18633-1' : -9.17,
                 'lig_18658-1' : -9.7,
                 'lig_18638-1' : -10.09,
                 'lig_18625-1' : -8.11,
                 'lig_18659-1' : -9.47,
                 'lig_18626-1' : -8.87}

    eps = 0.01
    for key, item in jacs_data.items():
        assert pytest.approx(item, eps) == df[df.name == key][('DerivedMeasurement', 'dg')].values[0].value_in_unit(unit.kilocalories_per_mole)
