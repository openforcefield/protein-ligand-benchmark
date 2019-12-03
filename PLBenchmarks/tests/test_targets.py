"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import ligands
import pytest
import yaml
import os
try:
    from importlib.resources import open_text, contents, is_resource
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text, contents, is_resource



def test_targets():
    file = open_text('PLBenchmarks.data', 'targets.yml')
    target_list = yaml.full_load(file)
    for target in target_list:
        # check if target directory is available
        assert target['dir'] in contents('PLBenchmarks.data') 
        assert not is_resource('PLBenchmarks.data', target['dir'])

        # check if YAML files of target are available
        assert 'target.yml' in contents(f'PLBenchmarks.data.{target["dir"]}.00_data') 
        assert 'ligands.yml' in contents(f'PLBenchmarks.data.{target["dir"]}.00_data') 
        
