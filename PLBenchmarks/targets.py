"""
targets.py
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
Target Data handling
"""

import yaml
from PLBenchmarks import ligands, edges

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

from PIL import Image

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
    for td in target_list:
        if td['name'] == target:
            return td['dir']
            break
    else:
        print('Directory for target {target} not found.')

def getTargetDataPath(target):
    print(target_list)
    for td in target_list:
        if td['name'] == target:
            return ['PLBenchmarks', 'data', td['dir'], '00_data']
            break
    else:
        print(f'Path for target {target} not found.')

class target:
    def __init__(self, name: str):
        """
            Initialize edge class from dictionary
        """
        self._name = name
        path = getTargetDataPath(self._name)
        file = open_text('.'.join(path), 'target.yml')
        data = yaml.full_load(file)
        self.data = pd.Series(data)
        self._ligands = None
        self._edges = None
        

    def getName(self):
        return self._name

    def getLigands(self):
        if self._ligands is None:
            self._ligands = ligands.getLigandSet(self._name)
        return self._ligands

    def getEdges(self):
        if self._edges is None:
            self._edges = edges.getEdgesDF(self._name)
        return self._edges
            
    def getDF(self, cols=None):
        if cols:
            return self.data[cols]
        else:
            return self.data

    def getGraph(self):
        ligs = self.getLigands()
        edgs = self.getEdges()
        
        G = nx.Graph()

        for i, item in enumerate(ligs):
            G.add_node(item.getName().split('_')[1], image=item.getImg())
        G.add_edges_from(edgs.values)

        pos=nx.circular_layout(G)

        fig = plt.figure(figsize=(40,20))
        ax = fig.gca()
        nx.draw(G,pos, node_size=35000, ax=ax, node_color=[[1,1,1,0]])

        label_pos = 0.5 # middle of edge, halfway between nodes
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


def getTargetsDF():
    dfs = []
    for td in target_list:
        dfs.append(target(td['name']).getDF())
    return pd.DataFrame(dfs)    
