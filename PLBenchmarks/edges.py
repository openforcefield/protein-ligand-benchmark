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
    Store and convert the data of one perturbation ("edge") in a `pandas.Series`.
    :param d: dict with the edge data
    :return None
    """
    def __init__(self, d: dict):
        """
        Initialize edge class from a dictionary
        :param d: dict with the edge data
        :return: None
        """
        self.data = pd.Series(d)
        self._name = None


    def addLigData(self, ligs):
        """
        Adds data from ligands to `edge`. Molecule images and the affinity difference are added.
        :param ligs: ligandSet class of the same target
        :return: None
        """
        for key, item in ligs.items():
            if key == 'lig_' + str(self.data[0]):
                l0 = item.data['ROMol'][0][0]
                dg0 = item.data[('DerivedMeasurement', 'dg')]
            if key == 'lig_' + str(self.data[1]):
                l1 = item.data['ROMol'][0][0]
                dg1 = item.data[('DerivedMeasurement', 'dg')]
        self.data['Mol1'] = l0
        self.data['Mol2'] = l1
        self.data['exp. DeltaG [kcal/mol]'] = (dg1-dg0).value_in_unit(unit.kilocalories_per_mole).round(2)

    def getDF(self, cols=None):
        """
        Access the edge data as a `pandas.DataFrame`
        :param cols: list of columns which should be returned in the `pandas.DataFrame`
        :return: data: `pandas.DataFrame`
        """
        if cols:
            return self.data[cols]
        else:
            return self.data

    def getDict(self):
        """
        Access the edge data as a dict which contains the name of the edge as key and the names of the two ligands as list.
        :return: edge: dict
        """
        return {f'edge_{self.data[0]}_{self.data[1]}': [f'lig_{self.data[0]}', f'lig_{self.data[1]}']}

    def getName(self):
        """
        Access the name of the edge.
        :return: name: string
        """
        if self._name is not None:
            return self._name
        else:
            return f'edge_{self.data[0]}_{self.data[1]}'
    

class edgeSet(dict):
    """
    Class inherited from dict to store all available edges of one target.
    """
    
    def __init__(self, target, *arg,**kw):
        """
        Initializes edgeSet class
        :param target: string name of target
        :param arg: arguments for dict (base class)
        :param kw: keywords for dict (base class)
        """
        super(edgeSet, self).__init__(*arg, **kw)
        tp = targets.getTargetDataPath(target)
        ligs = ligands.ligandSet(target)
        file = open_text('.'.join(tp), 'edges.yml')
        data = yaml.full_load_all(file)
        for d in data:
            e = edge(d)
            e.addLigData(ligs)
            self[e.getName()] = e

    def getEdge(self, name):
        """
        Accesses one edge of the edgeSet
        :param name: string name of the edge
        :return: edge: edge class
        """
        for key in self.keys():
            if key == name:
                return self[key]
                break
        else:
            raise ValueError(f'Edge {name} not part of set.')

    def getDF(self, columns=None):
        """
        Access the edgeSet as a `pandas.DataFrame`
        :param cols: list of columns which should be returned in the `pandas.DataFrame`
        :return: df: `pandas.DataFrame`
        """
        dfs=[]
        for key, item in self.items():
            dfs.append(item.getDF(columns))
        df = pd.DataFrame(dfs)
        return df

    def getHTML(self, columns=None):
        """
        Access the edgeSet as a HTML string
        :param cols: list of columns which should be returned in the `pandas.DataFrame`
        :return: html: HTML string
        """
        df = self.getDF(columns)
        html = df.to_html()
        return html

    def getDict(self):
        """
        Access the edgeSet as a dict which contains the name of the edges as key and the names of the two ligands in a list as items.
        :return: res: dict
        """
        res = {}
        for key, item in self.items():
            res.update(item.getDict())
        return res
        
