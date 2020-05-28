"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import targets
import pytest
import yaml
import os

try:
    from importlib.resources import open_text, contents, is_resource
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text, contents, is_resource


def test_targets():
    for target in targets.target_list:
        # check if target directory is available
        assert target["dir"] in contents("PLBenchmarks.data")
        assert not is_resource("PLBenchmarks.data", target["dir"])

        # check if YAML files of target are available
        assert "target.yml" in contents(f'PLBenchmarks.data.{target["dir"]}.00_data')
        assert "ligands.yml" in contents(f'PLBenchmarks.data.{target["dir"]}.00_data')


@pytest.mark.parametrize(
    "targetName, target", [(target["name"], target) for target in targets.target_list]
)
def test_target_class(targetName, target):
    tgt = targets.target(target["name"])
    tgt.getGraph()
    assert tgt.getName() == target["name"]


def test_targetSet():
    tgts = targets.targetSet()
    tgts.getDF()
    tgts.getHTML()
