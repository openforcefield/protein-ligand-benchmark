"""
targets.py
Functions and classes for handling the target data.
"""

import yaml
from PLBenchmarks import ligands, edges, utils

import numpy as np

import matplotlib.pyplot as plt

import pandas as pd

import networkx as nx

try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text


file = open_text('PLBenchmarks.data', 'targets.yml')
target_list = yaml.full_load(file)

def getTargetDir(target):
    """
    Gets the directory name of the target

    :param target: string with target name
    :return: string with  directory name
    """
    for td in target_list:
        if td['name'] == target:
            return td['dir']
            break
    else:
        print('Directory for target {target} not found.')

def getTargetDataPath(target):
    """
    Gets the file path of the target data

    :param target: string with target name
    :return: list of directories (have to be joined with '/' to get the file path relative to the PLBenchmarks repository)

    """
    for td in target_list:
        if td['name'] == target:
            return ['PLBenchmarks', 'data', td['dir'], '00_data']
            break
    else:
        raise ValueError(f'Path for target {target} not found.')

class target:
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
        path = getTargetDataPath(self._name)
        file = open_text('.'.join(path), 'target.yml')
        data = yaml.full_load(file)
        self.data = pd.Series(data)
        self.ligData = None
        self.htmlData = None
        self._ligands = None
        self._edges = None
        

    def getName(self):
        """
        Access the name of the target.

        :return: name as a string
        """
        return self._name

    def getLigandSet(self):
        """
        Get :py:class:`~PLBenchmarks.ligands.ligandSet` associated with the target

        :return: :py:class:`PLBenchmarks.ligands.ligandSet` object
        """
        if self._ligands is None:
            self._ligands = ligands.ligandSet(self._name)
        return self._ligands

    def addLigandData(self):
        """
        Adds data from ligands to :py:class:`PLBenchmarks.targets.target`. Molecule images and the minimum and maximum affinity are added.

        :return: None
        """
        lgs = self.getLigandSet()
        self.ligData = pd.Series({'numLigands': len(lgs)})
        affinities = []
        for key, item in lgs.items():
            affinities.append(item.data[('DerivedMeasurement', 'dg')].to_base_units().magnitude)
        self.ligData['maxDG'] = round(max(affinities), 1) * utils.ureg('kJ / mole')
        self.ligData['minDG'] = round(min(affinities), 1) * utils.ureg('kJ / mole')
        # calculation of the mean absolute deviation
        mean = np.average(affinities)
        mad = np.average(np.fabs(affinities - mean)) * utils.ureg('kJ / mole')
        self.ligData['MAD(DG)'] = round(mad, 1)

    def getLigData(self):
        if self.ligData is None:
            self.addLigandData()
        return self.ligData

    def getLigandSetDF(self, columns=None):
        """
        Get :py:class:`~PLBenchmarks.ligands.ligandSet` associated with the target in a :py:class:`pandas.DataFrame`

        :param columns: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        return self.getLigandSet().getDF(columns)

    def getLigandSetHTML(self, columns=None):
        """
        Get :py:class:`~PLBenchmarks.ligands.ligandSet` associated with the target in a html string

        :param columns: list of columns which should be returned
        :return: html string
        """
        return self.getLigandSet().getHTML(columns)
    
    def getEdgeSet(self):
        """
        Get :py:class:`~PLBenchmarks:edges:edgeSet` associated with the target

        :return: :py:class:`PLBenchmarks:edges:edgeSet` object
        """
        if self._edges is None:
            self._edges = edges.edgeSet(self._name)
        return self._edges

    def getEdgeSetDF(self, columns=None):
        """
        Get :py:class:`~PLBenchmarks:edges:edgeSet` associated with the target as a :py:class:`pandas.DataFrame`

        :param columns: list of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`PLBenchmarks:edges:edgeSet` object
        """
        return self.getEdgeSet().getDF(columns)

    def getEdgeSetHTML(self, columns=None):
        """
        Get :py:class:`~PLBenchmarks:edges:edgeSet` associated with the target in a html string

        :param columns: :py:class:`list` of edge which should be returned
        :return: html string
        """
        return self.getEdgeSet().getHTML(columns)        
    
    def getDF(self, columns=None):
        """
        Access the target data as a :py:class:`pandas.DataFrame`

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return:  :py:class:`pandas.DataFrame`
        """
        df = self.data
        df = df.append(self.getLigData())
        df = df.append(self.getHtmlData())
        if columns:
            return df[columns]
        else:
            return df

    def findLinks(self):
        """
        Processes primary data to have links in the html string of the target data

        :return: None
        """
        self.htmlData = pd.Series()
        if 'references' in list(self.data.index):
#            self.data.index = pd.MultiIndex.from_arrays([list(self.data.index), ['' for i in self.data.index]])
            refs = self.data['references']
            for key, item in refs.items():
                res = []
                if item is None:
                    continue
                for doi in item:
                    if str(doi) != 'nan':
                        res.append(utils.findDoiUrl(doi))
                self.htmlData[key] = (r'\n').join(res)
        if ('pdb') in list(self.data.index):
            pdb = self.data['pdb']
            self.htmlData['pdblinks'] = utils.findPdbUrl(' '.join(pdb.split(',')))

    def getHtmlData(self):
        if self.htmlData is None:
            self.findLinks()
        return self.htmlData

    def getGraph(self):
        """
        Get a graph representation of the ligand perturbations associated with the target in a :py:class:`matplotlib.figure`

        :return: :py:class:`matplotlib.figure`
        """
        
        G = nx.Graph()

        for key, item in self.getLigandSet().items():
            G.add_node(key.split('_')[1], image=item.getImg())
        G.add_edges_from([[item[0].split('_')[1], item[1].split('_')[1]] for key, item in self.getEdgeSet().getDict().items()])
        pos=nx.circular_layout(G)

        fig = plt.figure(figsize=(40,20))
        ax = fig.gca()
        nx.draw(G, pos, node_size=35000, ax=ax, node_color=[[1,1,1,0]])

        trans = ax.transData.transform
        trans2 = fig.transFigure.inverted().transform
        imsize = 0.15 # this is the image size
        for n in G.nodes():
            (x,y) = pos[n]
            xx,yy = trans((x,y)) # figure coordinates
            xa,ya = trans2((xx,yy)) # axes coordinates
            img =  G.nodes[n]['image']
            a = plt.axes([xa-imsize/2.0,ya-imsize/2.0, imsize, imsize ], fc=(1,1,1,0.0))
            a.set_xticks([]) 
            a.set_yticks([]) 
            a.spines['right'].set_visible(False)
            a.spines['top'].set_visible(False)
            a.spines['bottom'].set_visible(False)
            a.spines['left'].set_visible(False)
            a.imshow(img, alpha=1)
            a.set_aspect('equal')
        a.axis('off')
        return fig


class targetSet(dict):
    """
    Class inherited from dict to store all available targets in PLBenchmarks.
    """
    
    def __init__(self, *arg,**kw):
        """
        Initializes the :py:class:`~targets.targetSet` class

        :param target: string name of target
        :param arg: arguments for :py:class:`dict` (base class)
        :param kw: keywords for :py:class:`dict` (base class)
        """
        super(targetSet, self).__init__(*arg, **kw)
        for td in target_list:
            tgt = target(td['name'])
            self[tgt.getName()] = tgt
        self._df = None
          
    def getTarget(self, name):
        """
        Accesses one target of the targetSet

        :param name: string name of the target
        :return: :py:class:`PLBenchmarks.targets.target` class
        """
        tgt = None
        for key in self.keys():
            if key == name:
                tgt = self[key]
                break
        else:
            raise ValueError(f'Target {name} not part of set.')
        return tgt

    def getDF(self, columns=None):
        """
        Convert targetSet class to :py:class:`pandas.DataFrame`

        :param columns: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: :py:class:`pandas.DataFrame`
        """
        if self._df is None:
            dfs=[]
            for key in self.keys():
                self[key].addLigandData()
                self[key].findLinks()
                dfs.append(self[key].getDF())
            df = pd.DataFrame(dfs)
            self._df =df

        if columns is None:
            return self._df
        elif all(item in list(self._df.index) for item in columns):
            return self._df[columns]
        else:
            for item in columns:
                if item not in list(self._df.columns):
                    raise ValueError(f'Column {item} is not known and cannot be generated.')

    def getHTML(self, columns=None):
        """
        Access the :py:class:`~PLBenchmarks:targets:targetSet` as a HTML string

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.DataFrame`
        :return: HTML string
        """
        df = self.getDF(columns=columns)
        html = df.to_html()
        html = html.replace('REP1', '<a target="_blank" href="')
        html = html.replace('REP2', '">')
        html = html.replace('REP3', '</a>')
        html = html.replace("\\n","<br>")
        return html

    def getNames(self):
        """
        Get a list of available target names

        :return: :py:class:`list` of strings
        """
        return [key for key in self.keys()]


