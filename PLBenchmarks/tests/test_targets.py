"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
import pytest
import os

import PLBenchmarks
from PLBenchmarks import targets




def test_targets():
    for target in targets.target_list:
        # check if target directory is available
        assert target["dir"] in os.listdir(os.path.join(PLBenchmarks.__path__[0], "sample_data"))

        # check if YAML files of target are available
        assert "target.yml" in os.listdir(os.path.join(PLBenchmarks.__path__[0], "sample_data", target["dir"], "00_data"))
        assert "ligands.yml" in os.listdir(os.path.join(PLBenchmarks.__path__[0], "sample_data", target["dir"], "00_data"))


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
