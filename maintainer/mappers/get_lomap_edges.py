import pathlib
import yaml
import numpy as np
from rdkit import Chem
from openfe.setup import SmallMoleculeComponent
from openfe.setup.ligand_network_planning import minimal_spanning_graph
from openfe.setup.lomap_mapper import LomapAtomMapper
from openfe.setup.lomap_scorers import default_lomap_score


# Pick up ligand paths, load in rdkit and feed to SmallMoleculeComponent
path = pathlib.Path('../02_ligands/')

ligand_mols = []  # List of SmallMoleculeComponents

for p in path.rglob('*/crd/*sdf'):
    mol = [entry for entry in Chem.SDMolSupplier(str(p), removeHs=False)]
    ligand_mols.append(SmallMoleculeComponent(mol[0]))

# Create the minimum spanning graph

lomap_mapper = LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1)

network = minimal_spanning_graph(ligand_mols,
                                 [lomap_mapper,],
                                 default_lomap_score)

edges = [edge for edge in network.edges]

# Analyze distance overlap for manual validation
for i, edge in enumerate(edges):
    if np.any(edge.get_distances() > 1):
        name = f"{edge.molA.name}_{edge.molB.name}"
        print(f"{name} exceeds 1 A pairwise distance")

output = {}
output['mapper'] = "Lomap v2.1.0 MCS (time=20, threed=True, element_changes=False, max3d=1), default scorer"
output['planner'] = "openfe v0.4-dev minimal_spanning_graph"

edges_dict = {}

for edge in network.edges:
    edge_name = f'edge_{edge.molA.name}_{edge.molB.name}'
    mapping = list(edge.molA_to_molB.items())
    edges_dict[edge_name] = {
        'ligand_a': edge.molA.name,
        'ligand_b': edge.molB.name,
        'atom mapping': ' '.join(str(i) for i in mapping),
        'score': f"{edge.annotations['score']:.3f}",}

output['edges'] = edges_dict

with open('01_lomap_openfe.yaml', 'w') as yaml_file:
    yaml.dump(output, yaml_file, default_flow_style=False, sort_keys=False)
