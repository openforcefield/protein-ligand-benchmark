"""
targets.py
Functions and classes for handling the target data.
"""

import yaml
from PLBenchmarks import ligands, edges, utils

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
            Store and convert the data of one target in a :func:`pandas.Series`.
        :param name: string with target name
        :return: None
        """
        
        self._name = name
        path = getTargetDataPath(self._name)
        file = open_text('.'.join(path), 'target.yml')
        data = yaml.full_load(file)
        self.data = pd.Series(data)
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
         Get :func:`~PLBenchmarks.ligands.ligandSet` associated with the target
        :return: :func:`PLBenchmarks.ligands.ligandSet` object
        """
        if self._ligands is None:
            self._ligands = ligands.ligandSet(self._name)
        return self._ligands

    def addLigandData(self):
        """
         Adds data from ligands to :func:`PLBenchmarks.targets.target`. Molecule images and the minimum and maximum affinity are added.
        :return: None
        """
        lgs = self.getLigandSet()
        self.data['numLigands'] = len(lgs)
        affinities = []
        for key, item in lgs.items():
            affinities.append(item.data[('DerivedMeasurement', 'dg')])
        self.data['maxDG'] = max(affinities)
        self.data['minDG'] = min(affinities)        

    def getLigandSetDF(self, columns=None):
        """
         Get :func:`~PLBenchmarks.ligands.ligandSet` associated with the target in a :func:`pandas.DataFrame`
        :param columns: :func:`list` of columns which should be returned in the :func:`pandas.DataFrame`
        :return: :func:`pandas.DataFrame`
        """
        return self.getLigandSet().getDF(columns)

    def getLigandSetHTML(self, columns=None):
        """
         Get :func:`~PLBenchmarks.ligands.ligandSet` associated with the target in a html string
        :param columns: list of columns which should be returned
        :return: html string
        """
        return self.getLigandSet().getHTML(columns)
    
    def getEdgeSet(self):
        """
         Get :func:`~PLBenchmarks:edges:edgeSet` associated with the target
        :return: :func:`PLBenchmarks:edges:edgeSet` object
        """
        if self._edges is None:
            self._edges = edges.edgeSet(self._name)
        return self._edges

    def getEdgeSetDF(self, columns=None):
        """
         Get :func:`~PLBenchmarks:edges:edgeSet` associated with the target as a :func:`pandas.DataFrame`
        :param columns: list of columns which should be returned in the :func:`pandas.DataFrame`
        :return: :func:`PLBenchmarks:edges:edgeSet` object
        """
        return self.getEdgeSet().getDF(columns)

    def getEdgeSetHTML(self, columns=None):
        """
         Get :func:`~PLBenchmarks:edges:edgeSet` associated with the target in a html string
        :param columns: :func:`list` of edge which should be returned
        :return: html string
        """
        return self.getEdgeSet().getHTML(columns)        
    
    def getDF(self, columns=None):
        """
         Access the target data as a :func:`pandas.DataFrame`
        :param cols: :func:`list` of columns which should be returned in the :func:`pandas.DataFrame`
        :return:  :func:`pandas.DataFrame`
        """
        if columns:
            return self.data[columns]
        else:
            return self.data

    def findLinks(self):
        """
         Processes primary data to have links in the html string of the target data
        :return: None
        """
        if ('measurement', 'doi') in list(self.data.index):
            doi = self.data['measurement', 'doi']
            if str(doi) != 'nan':
                res = []
                for ddoi in re.split(r'[; ]+', str(doi)):
                    res.append(utils.findDoiUrl(ddoi))
            self.data['measurement', 'doi_html'] = (r'\n').join(res)
        if ('pdb') in list(self.data.index):
            pdb = self.data['pdb']
            self.data['pdb_html'] = utils.findPdbUrl(' '.join(pdb.split(',')))
        self.data.drop(columns=['pdb'], inplace=True)
        self.data.rename(columns={'pdb_html': 'pdb'}, inplace=True)
        
    def getGraph(self):
        """
         Get a graph representation of the ligand perturbations associated with the target in a :func:`matplotlib.figure`
        :return: :func:`matplotlib.figure`
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
         Initializes the :func:`~targets.targetSet` class
        :param target: string name of target
        :param arg: arguments for :func:`dict` (base class)
        :param kw: keywords for :func:`dict` (base class)
        """
        super(targetSet, self).__init__(*arg, **kw)
        for td in target_list:
            tgt = target(td['name'])
            tgt.findLinks()
            tgt.addLigandData()
            self[tgt.getName()] = tgt
          
    def getTarget(self, name):
        """
         Accesses one target of the targetSet
        :param name: string name of the target
        :return: :func:`PLBenchmarks.targets.target` class
        """
        for key in self.keys():
            if key == name:
                return self[key]
                break
        else:
            raise ValueError(f'Target {name} not part of set.')

    def getDF(self, columns=None):
        """
         Convert targetSet class to :func:`pandas.DataFrame`
        :param columns: :func:`list` of columns which should be returned in the :func:`pandas.DataFrame`
        :return: :func:`pandas.DataFrame`
        """
        dfs=[]
        for key, item in self.items():
            dfs.append(item.getDF(columns=columns))
        df = pd.DataFrame(dfs)
        return df

    def getHTML(self, columns=None):
        """
         Access the :func:`~PLBenchmarks:targets:targetSet` as a HTML string
        :param cols: :func:`list` of columns which should be returned in the :func:`pandas.DataFrame`
        :return: HTML string
        """
        df = self.getDF(columns)
        html = df.to_html()
        html = html.replace('REP1', '<a target="_blank" href="')
        html = html.replace('REP2', '">')
        html = html.replace('REP3', '</a>')
        html = html.replace("\\n","<br>")
        return html

    def getNames(self):
        """
         Get a list of available target names
        :return: :func:`list` of strings
        """
        return [key for key in self.keys()]


