"""
edges.py
Functions and classes for handling the perturbation edges.
"""
import os
import pandas as pd
import numpy as np
import yaml

from . import targets, ligands


class Edge:
    """
    Store and convert the data of one perturbation ("edge") in a :py:class:`pandas.Series`.

    :param d: :py:class:`dict` with the edge data
    :return: None
    """

    def __init__(self, d: dict):
        """
        Initialize edge class from a dictionary

        :param d: :py:class:`dict` with the edge data
        :return: None
        """
        self._data = pd.Series(d)
        self._name = None

    def add_ligand_data(self, ligand_set):
        """
        Adds data from ligands to :py:class:`~PLBenchmarks:edges.edge`. Molecule images and the affinity difference are added.

        :param ligand_set: :py:class:`PLBenchmarks:ligands:ligandSet` class of the same target
        :return: None
        """
        ligand0 = None
        ligand1 = None
        smiles0 = None
        smiles1 = None
        delta_g0 = 0
        delta_g1 = 0
        error0 = 0
        error1 = 0
        for key, item in ligand_set.items():
            if key == "lig_" + str(self._data[0]):
                ligand0 = item._data["ROMol"][0]
                delta_g0 = item._data[("DerivedMeasurement", "value")]
                smiles0 = item._data["smiles"]
                error0 = item._data[("DerivedMeasurement", "error")]
            if key == "lig_" + str(self._data[1]):
                ligand1 = item._data["ROMol"][0]
                smiles1 = item._data["smiles"]
                delta_g1 = item._data[("DerivedMeasurement", "value")]
                error1 = item._data[("DerivedMeasurement", "error")]
        self._data["Mol1"] = ligand0
        self._data["Mol2"] = ligand1
        self._data["Smiles1"] = smiles0
        self._data["Smiles2"] = smiles1
        self._data["exp. DeltaG [kcal/mol]"] = round(delta_g1 - delta_g0, 2)
        self._data["exp. Error [kcal/mol]"] = round(
            np.sqrt(np.power(error0, 2.0) + np.power(error1, 2.0)), 2
        )

    def get_dataframe(self, columns=None):
        """
        Access the edge data as a :py:class:`pandas.DataFrame`

        :param cols: list of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        if columns:
            return self._data[columns]
        else:
            return self._data

    def get_dict(self):
        """
        Access the edge data as a :py:class:`dict` which contains the name of the edge as key and the names of the two ligands as :py:class:`list`.

        :return: :py:class:`dict`
        """
        return {
            f"edge_{self._data[0]}_{self._data[1]}": [
                f"lig_{self._data[0]}",
                f"lig_{self._data[1]}",
            ]
        }

    def get_name(self):
        """
        Access the name of the edge.

        :return: name as string
        """
        if self._name is not None:
            return self._name
        else:
            return f"edge_{self._data[0]}_{self._data[1]}"


class EdgeSet(dict):
    """
    Class inherited from dict to store all available edges of one target.
    """

    def __init__(self, target, *arg, **kw):
        """
        Initializes edgeSet class

        :param target: string name of target
        :param arg: arguments for :py:class:`dict` (base class)
        :param kw: keywords for :py:class:`dict` (base class)
        """
        super(EdgeSet, self).__init__(*arg, **kw)
        target_path = targets.get_target_data_path(target)
        ligand_set = ligands.LigandSet(target)
        file = open(os.path.join(target_path, "edges.yml"))
        data = yaml.full_load_all(file)
        for d in data:
            edg = Edge(d)
            edg.add_ligand_data(ligand_set)
            self[edg.get_name()] = edg
        file.close()

    def get_edge(self, name):
        """
        Accesses one edge of the :py:class:`PLBenchmarks.edges.edgeSet`

        :param name: string name of the edge
        :return: :py:class:`PLBenchmarks:edges:edge` class
        """
        if name in self:
            return self[name]
        else:
            raise ValueError(f"Edge {name} not part of set.")

    def get_dataframe(self, columns=None):
        """
        Access the :py:class:`PLBenchmarks:edges.edgeSet` as a :py:class:`pandas.DataFrame`

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        dfs = []
        for key, item in self.items():
            dfs.append(item.get_dataframe(columns))
        df = pd.DataFrame(dfs)
        return df

    def get_html(self, columns=None):
        """
        Access the :py:class:`PLBenchmarks:edges.edgeSet` as a HTML string

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: HTML string
        """
        df = self.get_dataframe(columns)
        html_string = df.to_html()
        return html_string

    def get_dict(self):
        """
        Access the :py:class:`PLBenchmarks:edges.edgeSet` as a dict which contains the name of the edges as key and the names of the two ligands in a list as items.

        :return: :py:class:`dict`
        """
        result = {}
        for key, item in self.items():
            result.update(item.get_dict())
        return result
