"""
ProteinLigandBenchmarks
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
"""

# Add imports here
from .plbenchmarks import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
