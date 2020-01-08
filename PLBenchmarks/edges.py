"""
edges.py
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
Handles the perturbation edges
"""

from PLBenchmarks import utils, targets

import re
import pandas as pd
from simtk import unit
from rdkit.Chem import PandasTools
import yaml
try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text



class edge:
    
    """
        Store and convert the data of one perturbation ("edge") in a pandas Series.
    """
    
    def __init__(self, d: dict):
        """
            Initialize edge class from dictionary
        """
        self.data = pd.Series(d)


    def getDF(self, cols=None):
        if cols:
            return self.data[cols]
        else:
            return self.data

    def getDict(self):
        return {f'edge_{self.data[0]}_{self.data[1]}': [f'lig_{self.data[0]}', f'lig_{self.data[1]}']}
        
def getEdgesDF(target, cols=None):

    """
        Convenience function which returns all available ligands of one target in a `pandas` `dataframe`
    """

    tp = targets.getTargetDataPath(target)
    file = open_text('.'.join(tp), 'edges.yml')
    data = yaml.full_load_all(file)    
    dfs=[]
    for d in data:
        e = edge(d)
        dfs.append(e.getDF(cols))
    if cols:
        df = pd.DataFrame(dfs)[cols]
    else:
        df = pd.DataFrame(dfs)
    file.close()
    return df


def getEdgesDict( target ):
    tp = targets.getTargetDataPath(target)
    file = open_text('.'.join(tp), 'edges.yml')
    data = yaml.full_load_all(file)
    print(data)
    res = {}
    for d in data:
        e = edge(d)
        res.update(e.getDict())
    return res
