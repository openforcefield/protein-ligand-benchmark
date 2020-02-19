"""
ProteinLigandBenchmarks
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
"""

# Add imports here
import yaml

try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

file = open_text("PLBenchmarks.data", "targets.yml")
target_list = yaml.full_load(file)
