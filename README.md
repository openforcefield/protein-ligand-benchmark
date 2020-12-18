ProteinLigandBenchmarks
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/openforcefield/PLBenchmarks.svg?branch=master)](https://travis-ci.org/openforcefield/PLBenchmarks)
[![codecov](https://codecov.io/gh/openforcefield/PLBenchmarks/branch/master/graph/badge.svg)](https://codecov.io/gh/openforcefield/PLBenchmarks)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/openforcefield/PLBenchmarks.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/openforcefield/PLBenchmarks/context:python)
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
│   │   │   ├── protein.pdb                   #             aminoacid residues   
│   │   │   └── cofactors_crystalwater.pdb    #             cofactors and cyrstal waters    
│   │   └── top                               #         topology(s)
│   │   │   └── amber99sb-star-ildn-mut.ff    #             force field spec.           
│   │   │       ├── topol.itp                 #                 Gromacs ITP file
│   │   │       └── topol.top                 #                 Gromacs TOP file
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
- `date` is the data when the target was initially added to the benchmark set.
- `dir` is the directory name where all the data for the target is found. Usually it is the `date` and the `name` field, connected by a underscore `_`. 

#### `target.yml`

This file is always found in the meta data directory of each target: `<date>_<target_name>/00_data/target.yml`. It contains additionally information about the target:

#### `ligands.yml`

#### `edges.yml`

## Release History

Releases follow the `major.minor.micro` scheme recommended by [PEP440](https://www.python.org/dev/peps/pep-0440/#final-releases), where
- `major` increments denote a change that may break API compatibility with previous major releases and mark versions used in publications
- `minor` increments denote addition of new targets, changes of coordinates or topologies, or addition of new features to the API
- `micro` increments denote bugfixes and changes of metadata

## License

MIT. See the [License File](LICENSE) for more information.

CC-BY-4.0 for data (content of directory [`data`](data). See the [License File](LICENSE_DATA) for more information.

## Copyright

Copyright (c) 2019, Open Force Field Consortium, David F. Hahn


## Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
