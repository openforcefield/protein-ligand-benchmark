"""
edges.py
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
Handles the perturbation edges
"""

from PLBenchmarks import utils, targets, ligands

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


    def addLigData(self, ligs):
        for l in ligs:
            if l.data['name'][0] == 'lig_' + self.data[0]:
                l0 = l.data['ROMol'][0][0]
                dg0 = l.data[('DerivedMeasurement', 'dg')]
            if l.data['name'][0] == 'lig_' + self.data[1]:
                l1 = l.data['ROMol'][0][0]
                dg1 = l.data[('DerivedMeasurement', 'dg')]
        self.data['Mol1'] = l0
        self.data['Mol2'] = l1
        self.data['exp. DeltaG [kcal/mol]'] = (dg1-dg0).value_in_unit(unit.kilocalories_per_mole).round(2)

    def getDF(self, cols=None):
        if cols:
            return self.data[cols]
        else:
            return self.data

    def getDict(self):
        return {f'edge_{self.data[0]}_{self.data[1]}': [f'lig_{self.data[0]}', f'lig_{self.data[1]}']}


def getEdgesSet(target):

    """
        Convenience function which returns all available edges of one target in a list
    """

    tp = targets.getTargetDataPath(target)
    ligs = ligands.getLigandSet(target)    
    file = open_text('.'.join(tp), 'edges.yml')
    data = yaml.full_load_all(file)
    edgs = []
    for d in data:
        e = edge(d)
        e.addLigData(ligs)
        edgs.append(e)
    file.close()
    return edgs

def getEdgesSetDF(target, cols=None):

    """
        Convenience function which returns all available ligands of one target in a `pandas` `dataframe`
    """

    tp = targets.getTargetDataPath(target)
    ligs = ligands.getLigandSet(target)
    file = open_text('.'.join(tp), 'edges.yml')
    data = yaml.full_load_all(file)    
    dfs=[]
    for d in data:
        e = edge(d)
        e.addLigData(ligs)
        dfs.append(e.getDF(cols))
    if cols:
        df = pd.DataFrame(dfs)[cols]
    else:
        df = pd.DataFrame(dfs)
    file.close()
    return df


def getEdgesSetHTML(target, cols=None):

    """
        Convenience function which returns all edges of one target in an html string
    """

    df = getEdgesSetDF(target, cols=cols)

    html = df.to_html()
    return html

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
