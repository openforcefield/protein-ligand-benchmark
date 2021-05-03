ProteinLigandBenchmarks
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/openforcefield/protein-ligand-benchmark.svg?branch=master)](https://travis-ci.org/openforcefield/protein-ligand-benchmark)
[![codecov](https://codecov.io/gh/openforcefield/protein-ligand-benchmark/branch/master/graph/badge.svg)](https://codecov.io/gh/openforcefield/protein-ligand-benchmark)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/openforcefield/protein-ligand-benchmark.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/openforcefield/protein-ligand-benchmark/context:python)
[![Documentation Status](https://readthedocs.org/projects/plbenchmarks/badge/?version=latest)](https://plbenchmarks.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.

## Documentation

[Documentation](https://plbenchmarks.readthedocs.io/en/latest/) for the `openforcefield` toolkit is hosted at [readthedocs](https://plbenchmarks.readthedocs.io/en/latest/).

## Data file tree and file description

The data is organized as followed:  

```
data
├── targets.yml                               # list of all targets and their directories   
├── <date>_<target_name_1>                    # directory for target 1
│   ├── 00_data                               #     metadata for target 1
│   │   ├── edges.yml                         #         edges/perturbations
│   │   ├── ligands.yml                       #         ligands and activities
│   │   └── target.yml                        #         target
│   ├── 01_protein                            #     protein data
│   │   ├── crd                               #         coordinates
│   │   │   ├── cofactors_crystalwater.pdb    #             cofactors and cyrstal waters    
│   │   │   └── protein.pdb                   #             aminoacid residues   
│   │   └── top                               #         topology(s)
│   │   │   └── amber99sb-star-ildn-mut.ff    #             force field spec.     
│   │   │       ├── cofactors_crystalwater.top#                 Gromacs TOP file of cofactors and crystal water
│   │   │       ├── protein.top               #                 Gromacs TOP file of amino acid residues
│   │   │       └── *.itp                     #                 Gromacs ITP file(s) to be included in TOP files
│   └── 02_ligands                            #     ligands
│   ├── lig_<name_1>                          #          ligand 1 
│   │   ├── crd                               #              coordinates
│   │   │   └── lig_<name_1>.sdf              #                  SDF file
│   │   └── top                               #              topology(s)
│   │       └── openff-1.0.0.offxml           #                  force field spec.       
│   │           ├── fflig_<name_1>.itp        #                      Gromacs ITP file : atom types     
│   │           ├── lig_<name_1>.itp          #                      Gromacs ITP file       
│   │           ├── lig_<name_1>.top          #                      Gromacs TOP file                
│   │           └── posre_lig_<name_1>.itp    #                      Gromacs ITP file : position restraint file  
│   ├── lig_<name_2>                          #         ligand 2                               
│   …                                        
│   └── 03_hybrid                             #    edges (perturbations)
│   ├── edge_<name_1>_<name_2>                #         edge between ligand 1 and ligand 2   
│   │   └── water                             #             edge in water 
│   │       ├── crd                           #                 coordinates 
│   │       │   ├── mergedA.pdb               #                     merged conf based on coords of ligand 1  
│   │       │   ├── mergedB.pdb               #                     merged conf based on coords of ligand 2   
│   │       │   ├── pairs.dat                 #                     atom mapping                  
│   │       │   └── score.dat                 #                     similarity score         
│   │       └── top                           #                 topology(s)       
│   │           └── openff-1.0.0.offxml       #                     force field spec.         
│   │               ├── ffmerged.itp          #                         Gromacs ITP file  
│   │               ├── ffMOL.itp             #                         Gromacs ITP file   
│   │               └── merged.itp            #                         Gromacs ITP file     
│   …                                        
├── <date>_<target_name_2>                    # directory for target 2  
…
```
## Description of meta data YAML files

#### `targets.yml`

This file lists all the registered targets in the benchmark set. Each entry denotes one target and contains the following information:

```
mcl1_sample:
  name:     mcl1_sample
  date:     2020-08-26
  dir:      2020-08-26_mcl1_sample
```

`mcl1_sample` is the entry name and each entry has three sub-entries: 
- `name` is the target name, which is usually the same as the entry name of the target. 
- `date` is the date when the target was initially added to the benchmark set.
- `dir` is the directory name where all the data for the target is found. Usually it is the `date` and the `name` field, connected by a underscore `_`. 

#### `target.yml`

This file is found in the meta data directory of each target: `<date>_<target_name>/00_data/target.yml`. It contains additionally information about the target:

```
alternate:
  iridium_classifier: HT
  iridium_score: 0.3
  pdb: 6O6F
associated_sets:
- Schrodinger JACS
comments: hydrophobic interactions contributing to binding
date: 2019-12-13
dpi: 0.26
id: 9
iridium_classifier: HT
iridium_score: 0.41
name: mcl1
netcharge: 4 e
pdb: 4HW3
references:
  calculation:
  - 10.1021/ja512751q
  - 10.1021/acs.jcim.9b00105
  - 10.1039/C9SC03754C
  measurement:
  - 10.1021/jm301448p
```

Explanation of the entries:

- `alternate`: Alternate X-ray structure which could be used
  - `iridium_classifier`: Iridium classifier of the alternate structure
  - `iridium_score`: Iridium score of the alternate structure
  - `pdb`: PDB ID of the alternate structure
- `associated_sets`: list of benchmark set tags, where this target is in (e.g. `"Schrodinger JACS"`)
- `comments`: hydrophobic interactions contributing to binding
- `date`: date when the target was initially added to the benchmark set.
- `dpi`: diffraction precision index of the used structure (quality metric for the structure)
- `id`: a given ID
- `iridium_classifier`: Iridium classifier of the used structure
- `iridium_score`: Iridium score of the used structure
- `name`: name/identifier of the target
- `netcharge`: total charge of the prepared protein (this should be equalized with counter ions during preparation of the simulation system)
- `pdb`: PDB ID of the used structure
- `references`: doi to references
  - `calculation`: list of references where this target was used in calculations
  - `measurement`: list of references of affinity measurements
  
#### `ligands.yml`

This file is found in the meta data directory of each target: `<date>_<target_name>/00_data/ligands.yml`. It contains information of the ligands of one target. One entry looks like this:

```
lig_23:
  measurement:
    comment: Table 2, entry 23
    doi: 10.1021/jm301448p
    error: 0.03
    type: ki
    unit: uM
    value: 0.37
  name: lig_23
  smiles: '[H]c1c(c(c2c(c1[H])c(c(c(c2OC([H])([H])C([H])([H])C([H])([H])C3=C(Sc4c3c(c(c(c4[H])[H])[H])[H])C(=O)[O-])[H])[H])[H])[H])[H]'
```

Explanation of the entries:

- `measurement`: affinity measurement entry
  - `comment`: comment about the measurement
  - `doi`: DOI (digital object identifier) pointing to the reference for this measurement
  - `error`: Error of measurement, `null` if not reported
  - `type`: type of measurement observable, `ki` (binding equilibrium constant), `ic50` (IC50 value), `pic50` (pIC50 value), or `dg` (free energy of binding) are accepted entries. 
  - `unit`: Unit of value and error entries.
  - `value`: Value of the measurement.
- `name`: name of ligand, which always starts with `lig_`, followed by a unique identifier.
- `smiles`: SMILES string of the ligand, with charge state information and chirality information. 

#### `edges.yml`

This file is found in the meta data directory of each target: `<date>_<target_name>/00_data/edges.yml`. It contains information of the edges of one target.  One entry looks like this:

```
edge_50_60:
  ligand_a: lig_50
  ligand_b: lig_60
```

Each entry is just a list of two ligand identifiers. 

## Release History

Releases follow the `major.minor.micro` scheme recommended by [PEP440](https://www.python.org/dev/peps/pep-0440/#final-releases), where
- `major` increments denote a change that may break API compatibility with previous major releases
- `minor` increments denote addition of new targets or addition and larger changes to the API
- `micro` increments denote bugfixes, addition of API features, changes of coordinates or topologies, and changes of metadata

## License

MIT. See the [License File](LICENSE) for more information.

CC-BY-4.0 for data (content of directory [`data`](data). See the [License File](LICENSE_DATA) for more information.

## Contributions

- **Authors** David Hahn
- **Data Contributors** The authors of the following publications, expecially Vytautas Gapsys and Christina E. M. Schindler.
  - [V. Gapsys et al., Large scale relative protein ligand binding affinities using non-equilibrium alchemy, Chem. Sci., 2020,11, 1140-1152](https://doi.org/10.1039/C9SC03754C)
  - [Christina E. M. Schindler et al., Large-Scale Assessment of Binding Free Energy Calculations in Active Drug Discovery Projects, J. Chem. Inf. Model. 2020, 60, 11, 5457–5474](https://doi.org/10.1021/acs.jcim.0c00900)
  - [Laura Perez Benito et al., Predicting Activity Cliffs with Free-Energy Perturbation, J. Chem. Theory Comput. 2019, 15, 3, 1884–1895](https://pubs.acs.org/doi/10.1021/acs.jctc.8b01290)
- **Discussions and Suggestions** Christopher I. Bayly, Marko Breznik, Hannah E. Bruce Macdonald, John D.Chodera, Katharina Meier, Antonia S. J. S. Mey, David L. Mobley, Laura Perez Benito, Gary Tresadern, Gregory L. Warren and all members of the Open Force Field Initiative
- **Code review and discussions** Matt Thompson, Jeffrey Wagner

## Copyright

Copyright (c) 2019, Open Force Field Consortium, David F. Hahn

## Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
