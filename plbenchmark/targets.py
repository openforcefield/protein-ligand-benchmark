"""
targets.py
Functions and classes for handling the target data.
"""

import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

from . import __path__, ligands, edges, utils


data_path = os.path.abspath(os.path.join(os.path.join(__path__[0], "sample_data")))
file = open(os.path.join(data_path, "targets.yml"))
target_dict = yaml.full_load(file)
file.close()


def set_data_dir(path=os.path.abspath(os.path.join(__path__[0], "sample_data"))):
    """
    Gets the directory name of the target

    :param path: string with path to data directory
    """
    global data_path
    data_path = os.path.abspath(path)
    file = open(os.path.join(data_path, "targets.yml"))
    global target_dict
    target_dict = yaml.full_load(file)
    file.close()


def get_target_dir(target):
    """
    Gets the directory name of the target

    :param target: string with target name
    :return: string with  directory name
    """
    if target in target_dict:
        return target_dict[target]["dir"]
    else:
        raise ValueError(f"Directory for target {target} not found.")


def get_target_data_path(target):
    """
    Gets the file path of the target data

    :param target: string with target name
    :return: list of directories (have to be joined with '/' to get the file path relative to the plbenchmark repository)

    """
    if target in target_dict:
        return os.path.join(data_path, target_dict[target]["dir"], "00_data", "")
    else:
        raise ValueError(f"Path for target {target} not found.")


class Target:
    """
    Class to store the data of one target.

    """

    def __init__(self, name: str):
        """
        Store and convert the data of one target in a :py:class:`pandas.Series`.

        :param name: string with target name
        :return: None
        """

        self._name = name
        path = get_target_data_path(self._name)
        file = open(os.path.join(path, "target.yml"))
        data = yaml.full_load(file)
        self._data = pd.Series(data)
        file.close()

        self.ligand_data = None
        self.html_data = None
        self._ligands = None
        self._edges = None

    def get_name(self):
        """
        Access the name of the target.

        :return: name as a string
        """
        return self._name

    def get_ligand_set(self):
        """
        Get :py:class:`~plbenchmark.ligands.ligandSet` associated with the target

        :return: :py:class:`plbenchmark.ligands.ligandSet` object
        """
        if self._ligands is None:
            self._ligands = ligands.LigandSet(self._name)
        return self._ligands

    def add_ligand_data(self):
        """
        Adds data from ligands to :py:class:`plbenchmark.targets.target`. Molecule images and the minimum and maximum affinity are added.

        :return: None
        """
        lgs = self.get_ligand_set()
        self.ligand_data = pd.Series({"numLigands": len(lgs)})
        affinities = []
        for key, item in lgs.items():
            affinities.append(
                item._data[("DerivedMeasurement", "value")].to("kcal/mole").magnitude
            )
        self.ligand_data["maxDG"] = round(
            max(affinities) * utils.unit("kcal / mole"), 1
        )
        self.ligand_data["minDG"] = round(
            min(affinities) * utils.unit("kcal / mole"), 1
        )
        # calculation of the standard deviation
        std = np.std(affinities)
        self.ligand_data["std(DG)"] = round(std * utils.unit("kcal / mole"), 1)

    def get_ligand_data(self):
        if self.ligand_data is None:
            self.add_ligand_data()
        return self.ligand_data

    def get_ligand_set_dataframe(self, columns=None):
        """
        Get :py:class:`~plbenchmark.ligands.ligandSet` associated with the target in a :py:class:`pandas.DataFrame`

        :param columns: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        return self.get_ligand_set().get_dataframe(columns)

    def get_ligand_set_html(self, columns=None):
        """
        Get :py:class:`~plbenchmark.ligands.ligandSet` associated with the target in a html string

        :param columns: list of columns which should be returned
        :return: html string
        """
        return self.get_ligand_set().get_html(columns)

    def get_edge_set(self):
        """
        Get :py:class:`~plbenchmark:edges:edgeSet` associated with the target

        :return: :py:class:`plbenchmark:edges:edgeSet` object
        """
        if self._edges is None:
            self._edges = edges.EdgeSet(self._name)
        return self._edges

    def get_edge_set_dataframe(self, columns=None):
        """
        Get :py:class:`~plbenchmark:edges:edgeSet` associated with the target as a :py:class:`pandas.DataFrame`

        :param columns: list of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`plbenchmark:edges:edgeSet` object
        """
        return self.get_edge_set().get_dataframe(columns)

    def get_edge_set_html(self, columns=None):
        """
        Get :py:class:`~plbenchmark:edges:edgeSet` associated with the target in a html string

        :param columns: :py:class:`list` of edge which should be returned
        :return: html string
        """
        return self.get_edge_set().get_html(columns)

    def get_dataframe(self, columns=None):
        """
        Access the target data as a :py:class:`pandas.DataFrame`

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return:  :py:class:`pandas.DataFrame`
        """
        df = self._data
        df = df.append(self.get_ligand_data())
        df = df.append(self.get_html_data())
        if columns:
            return df[columns]
        else:
            return df

    def find_links(self):
        """
        Processes primary data to have links in the html string of the target data

        :return: None
        """
        self.html_data = pd.Series(dtype=object)
        if "references" in list(self._data.index):
            #            self._data.index = pd.MultiIndex.from_arrays([list(self._data.index), ['' for i in self._data.index]])
            refs = self._data["references"]
            for key, item in refs.items():
                res = []
                if item is None:
                    continue
                for doi in item:
                    if str(doi) != "nan":
                        res.append(utils.find_doi_url(doi))
                self.html_data[key] = (r"\n").join(res)
        if ("pdb") in list(self._data.index):
            pdb = self._data["pdb"]
            if pdb is None:
                self.html_data["pdblinks"] = ""
            else:
                self.html_data["pdblinks"] = utils.find_pdb_url(
                    " ".join(pdb.split(","))
                )

    def get_html_data(self):
        if self.html_data is None:
            self.find_links()
        return self.html_data

    def get_graph(self):
        """
        Get a graph representation of the ligand perturbations associated with the target in a :py:class:`matplotlib.figure`

        :return: :py:class:`matplotlib.figure`
        """

        graph = nx.Graph()

        for key, item in self.get_ligand_set().items():
            graph.add_node(key.split("_")[1], image=item.get_image())
        graph.add_edges_from(
            [
                [item["ligand_a"].split("_")[1], item["ligand_b"].split("_")[1]]
                for key, item in self.get_edge_set().get_dict().items()
            ]
        )
        pos = nx.circular_layout(graph)

        fig = plt.figure(figsize=(60, 40))
        ax = fig.gca()
        nx.draw(graph, pos, node_size=35000, ax=ax, node_color=[[1, 1, 1, 0]])

        trans = ax.transData.transform
        trans2 = fig.transFigure.inverted().transform
        imsize = 0.075  # this is the image size
        for n in graph.nodes():
            (x, y) = pos[n]
            xx, yy = trans((x, y))  # figure coordinates
            xa, ya = trans2((xx, yy))  # axes coordinates
            img = graph.nodes[n]["image"]
            a = plt.axes(
                [xa - imsize / 2.0, ya - imsize / 2.0, imsize, imsize],
                fc=(1, 1, 1, 0.0),
            )
            a.set_xticks([])
            a.set_yticks([])
            a.spines["right"].set_visible(False)
            a.spines["top"].set_visible(False)
            a.spines["bottom"].set_visible(False)
            a.spines["left"].set_visible(False)
            a.imshow(img, alpha=1)
            a.set_aspect("equal")
        a.axis("off")
        return fig


class TargetSet(dict):
    """
    Class inherited from dict to store all available targets in plbenchmark.
    """

    def __init__(self, *arg, **kw):
        """
        Initializes the :py:class:`~targets.targetSet` class

        :param target: string name of target
        :param arg: arguments for :py:class:`dict` (base class)
        :param kw: keywords for :py:class:`dict` (base class)
        """
        super(TargetSet, self).__init__(*arg, **kw)
        for name in target_dict.keys():
            target = Target(name)
            self[target.get_name()] = target
        self._df = None

    def __eq__(self, other):
        if not isinstance(other, TargetSet):
            return False
        return dict.__eq__(self, other) and self._df == other._df

    def __ne__(self, other):
        if not isinstance(other, TargetSet):
            return True
        return dict.__ne__(self, other) or self._df != other._df

    def get_target(self, name):
        """
        Accesses one target of the targetSet

        :param name: string name of the target
        :return: :py:class:`plbenchmark.targets.target` class
        """
        if name in self:
            return self[name]
        else:
            raise ValueError(f"Target {name} not part of set.")

    def get_dataframe(self, columns=None):
        """
        Convert targetSet class to :py:class:`pandas.DataFrame`

        :param columns: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        if self._df is None:
            dfs = []
            for key in self.keys():
                self[key].add_ligand_data()
                self[key].find_links()
                dfs.append(self[key].get_dataframe())
            df = pd.DataFrame(dfs)
            self._df = df

        if columns is None:
            return self._df
        elif all(item in list(self._df.columns) for item in columns):
            return self._df[columns]
        else:
            for item in columns:
                if item not in list(self._df.columns):
                    raise ValueError(
                        f"Column {item} is not known and cannot be generated."
                    )

    def get_html(self, columns=None):
        """
        Access the :py:class:`~plbenchmark:targets:targetSet` as a HTML string

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: HTML string
        """
        df = self.get_dataframe(columns=columns)
        html_string = df.to_html()
        html_string = html_string.replace("REP1", '<a target="_blank" href="')
        html_string = html_string.replace("REP2", '">')
        html_string = html_string.replace("REP3", "</a>")
        html_string = html_string.replace("\\n", "<br>")
        return html_string

    def get_names(self):
        """
        Get a list of available target names

        :return: :py:class:`list` of strings
        """
        return [key for key in self.keys()]
