"""
targets.py
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
Target Data handling
"""

import yaml
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
