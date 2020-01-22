Installing the Protein Ligand Benchmark Set
===========================================

The Protein Ligand Benchmark Set is currently only installable from source.

Installation from Source
------------------------

.. toctree::
   :maxdepth: 1

To install PLBenchmarks from source, clone the repository from `github
<https://github.com/openforcefield/PLBenchmarks>`_::

    git clone https://github.com/openforcefield/PLBenchmarks.git
    cd PLBenchmarks

Create a custom conda environment which contains the required dependencies and activate it::

    conda env create --name plbenchmarks --file devtools/conda-envs/test_env.yaml
    conda activate plbenchmarks

The final step is to install PLBenchmarks itself::

    python setup.py develop

