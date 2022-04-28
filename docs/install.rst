Installing the Protein Ligand Benchmark Set
===========================================

The Protein Ligand Benchmark Set is currently only installable from source.

Installation from Source
------------------------

.. toctree::
   :maxdepth: 1

The repository uses `git-lfs (large file storage) <https://git-lfs.github.com>`_ for the storage of all the data file. Ideally git-lfs is installed first before cloning the repository.

:: 

    conda create -n plbenchmark python=3.7 git-lfs
    conda activate plbenchmark
    git lfs clone https://github.com/openforcefield/protein-ligand-benchmark.git
    cd protein-ligand-benchmark
    conda env update --file environment.yml
    pip install -e .

Related Publication
-------------------
The `preprint <https://arxiv.org/pdf/2105.06222.pdf>`_ on "Best practices for constructing, preparing, and evaluating protein-ligand binding affinity benchmarks" provides accompanying information to this benchmark dataset and how to use it for alchemical free energy calculations.
For any suggestions of improvements please raise an issue in its `GitHub repository <https://github.com/openforcefield/FE-Benchmarks-Best-Practices>`_.
