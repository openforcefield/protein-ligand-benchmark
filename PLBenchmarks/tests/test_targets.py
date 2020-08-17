"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
import pytest
import os
import datetime
import pandas as pd
import PLBenchmarks
from PLBenchmarks import targets, ligands, edges


def test_targets():
    assert len(targets.target_list) == 1
    mcl1_dict = {
        "date": datetime.date(2019, 12, 13),
        "dir": "2019-12-13_mcl1",
        "name": "mcl1",
    }
    for key in mcl1_dict.keys():
        assert targets.target_list[0][key] == mcl1_dict[key]

    for target in targets.target_list:
        # check if target directory is available
        assert target["dir"] in os.listdir(
            os.path.join(PLBenchmarks.__path__[0], "sample_data")
        )
        assert target["dir"] == targets.getTargetDir(target["name"])

        # check if YAML files of target are available
        assert "target.yml" in os.listdir(
            os.path.join(
                PLBenchmarks.__path__[0], "sample_data", target["dir"], "00_data"
            )
        )
        assert "ligands.yml" in os.listdir(
            os.path.join(
                PLBenchmarks.__path__[0], "sample_data", target["dir"], "00_data"
            )
        )

    test_target = "aaa"
    with pytest.raises(ValueError, match=f"Path for target {test_target} not found."):
        targets.getTargetDataPath(test_target)
    test_target = "aaa"
    with pytest.raises(
        ValueError, match=f"Directory for target {test_target} not found."
    ):
        targets.getTargetDir(test_target)


def test_target_class():
    target = [t for t in targets.target_list if t["name"] == "mcl1"][0]
    tgt = targets.target(target["name"])
    assert tgt.getName() == target["name"]
    assert tgt.ligData == None
    assert tgt.htmlData == None
    assert tgt._ligands == None
    assert tgt._edges == None

    ligSet = ligands.ligandSet("mcl1")
    assert tgt.getLigandSet().keys() == ligSet.keys()
    assert tgt._ligands != None
    tgt.addLigandData()
    assert type(tgt.ligData) == type(pd.Series(dtype=object))
    pd.testing.assert_series_equal(tgt.ligData, tgt.getLigData())
    tgt.ligData = None
    ligData = tgt.getLigData()
    assert type(tgt.ligData) == type(pd.Series(dtype=object))
    assert ligData["numLigands"] == 42
    # cannot compare ROMol column (SVG image), that's why we only compare these columns
    columns = ["name", "smiles", "docked", "measurement", "DerivedMeasurement"]
    df1 = tgt.getLigandSetDF(columns=columns)
    df2 = ligSet.getDF(columns=columns)
    pd.testing.assert_frame_equal(df1, df2)
    assert tgt.getLigandSetHTML() == ligSet.getHTML()

    edgeSet = edges.edgeSet("mcl1")
    columns = [
        0,
        1,
        "Smiles1",
        "Smiles2",
        "exp. DeltaG [kcal/mol]",
        "exp. Error [kcal/mol]",
    ]
    pd.testing.assert_frame_equal(
        tgt.getEdgeSet().getDF(columns=columns), edgeSet.getDF(columns=columns)
    )
    assert tgt.getEdgeSet().getDict() == edgeSet.getDict()
    assert tgt.getEdgeSetHTML() == edgeSet.getHTML()

    # TODO: this actually does not test anything, only checks if it works
    tgt.findLinks()
    tgt.getHtmlData()
    tgt.getGraph()


def test_targetSet():
    tgts = targets.targetSet()
    tgts2 = targets.targetSet()
    assert tgts == tgts

    tgts2._df = pd.DataFrame()
    assert tgts != tgts2

    # TODO: implement __eq__ and __ne__ operator for target class
    # assert targets.target('mcl1') == tgts.getTarget('mcl1')

    df = tgts.getDF()
    assert df.shape[0] == 1
    assert df.shape[1] > 1
    assert df["name"][0] == "mcl1"

    df = tgts.getDF(columns=["name"])
    assert df.shape[0] == 1
    assert df.shape[1] == 1
    assert df["name"][0] == "mcl1"

    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        tgts.getDF(columns=["xxx"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        tgts.getDF(columns=["name", "xxx"])

    # TODO: this actually does not test anything, only checks if it works
    tgts.getHTML()
    tgts.getHTML(columns=["name"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        tgts.getHTML(columns=["xxx"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        tgts.getHTML(columns=["name", "xxx"])
