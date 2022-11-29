"""
Unit and regression test for the plbenchmark package.
"""

# Import package, test suite, and other packages as needed
import pytest
import os
import datetime
import pandas as pd
import plbenchmark
from plbenchmark import targets, ligands, edges


def test_targets():
    assert len(targets.target_dict) == 1
    mcl1_dict = {
        "mcl1_sample": {
            "date": datetime.date(2020, 8, 26),
            "dir": "2020-08-26_mcl1_sample",
            "name": "mcl1_sample",
        }
    }
    for key in mcl1_dict["mcl1_sample"].keys():
        assert targets.target_dict["mcl1_sample"][key] == mcl1_dict["mcl1_sample"][key]

    for target, item in targets.target_dict.items():
        # check if target directory is available
        assert item["dir"] in os.listdir(
            os.path.join(plbenchmark.__path__[0], "sample_data")
        )
        assert item["dir"] == targets.get_target_dir(target)

        # check if YAML files of target are available
        assert "target.yml" in os.listdir(
            os.path.join(
                plbenchmark.__path__[0], "sample_data", item["dir"], "00_data"
            )
        )
        assert "ligands.yml" in os.listdir(
            os.path.join(
                plbenchmark.__path__[0], "sample_data", item["dir"], "00_data"
            )
        )

    test_target = "aaa"
    with pytest.raises(ValueError, match=f"Path for target {test_target} not found."):
        targets.get_target_data_path(test_target)
    test_target = "aaa"
    with pytest.raises(
        ValueError, match=f"Directory for target {test_target} not found."
    ):
        targets.get_target_dir(test_target)


def test_target_class():
    target = targets.target_dict["mcl1_sample"]
    tgt = targets.Target(target["name"])
    assert tgt.get_name() == target["name"]
    assert tgt.ligand_data == None
    assert tgt.html_data == None
    assert tgt._ligands == None
    assert tgt._edges == None

    ligand_set = ligands.LigandSet("mcl1_sample")
    assert tgt.get_ligand_set().keys() == ligand_set.keys()
    assert tgt._ligands != None
    tgt.add_ligand_data()
    assert type(tgt.ligand_data) == type(pd.Series(dtype=object))
    pd.testing.assert_series_equal(tgt.ligand_data, tgt.get_ligand_data())
    tgt.ligand_data = None
    ligand_data = tgt.get_ligand_data()
    assert type(tgt.ligand_data) == type(pd.Series(dtype=object))
    assert ligand_data["numLigands"] == 15
    # cannot compare ROMol column (SVG image), that's why we only compare these columns
    columns = ["name", "smiles", "measurement", "DerivedMeasurement"]
    df1 = tgt.get_ligand_set_dataframe(columns=columns)
    df2 = ligand_set.get_dataframe(columns=columns)
    pd.testing.assert_frame_equal(df1, df2)
    # Temporarily disable in #82 - RDKit mol hash to different values
    #assert tgt.get_ligand_set_html() == ligand_set.get_html()

    edge_set = edges.EdgeSet("mcl1_sample")
    columns = [
        "ligand_a",
        "ligand_b",
        "Smiles1",
        "Smiles2",
        "exp. DeltaG [kcal/mol]",
        "exp. Error [kcal/mol]",
    ]
    pd.testing.assert_frame_equal(
        tgt.get_edge_set().get_dataframe(columns=columns),
        edge_set.get_dataframe(columns=columns),
    )
    dict1 = tgt.get_edge_set().get_dict()
    dict2 = edge_set.get_dict()
    assert dict1.keys() == dict2.keys()
    for key, item in dict1.items():
        assert item.keys() == dict2[key].keys()
        for kk, ii in item.items():
            if kk != "Mol1" and kk != "Mol2":
                assert ii == dict2[key][kk]
    # Temporarily disable in #82 - RDKit mol hash to different values
    #assert tgt.get_edge_set_html() == edge_set.get_html()

    # TODO: this actually does not test anything, only checks if it works
    tgt.find_links()
    tgt.get_html_data()
    tgt.get_graph()


def test_target_set():
    target_set = targets.TargetSet()
    target_set_2 = targets.TargetSet()
    # TODO: implement __eq__ and __ne__ operator for TargetSet class
    # assert target_set == target_set_2

    target_set_2._df = pd.DataFrame()
    assert target_set != target_set_2

    # TODO: implement __eq__ and __ne__ operator for Target class
    # assert targets.Target('mcl1_sample') == target_set.get_target('mcl1_sample')

    df = target_set.get_dataframe()
    assert df.shape[0] == 1
    assert df.shape[1] > 1
    assert df["name"][0] == "mcl1_sample"

    df = target_set.get_dataframe(columns=["name"])
    assert df.shape[0] == 1
    assert df.shape[1] == 1
    assert df["name"][0] == "mcl1_sample"

    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        target_set.get_dataframe(columns=["xxx"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        target_set.get_dataframe(columns=["name", "xxx"])

    # TODO: this actually does not test anything, only checks if it works
    target_set.get_html()
    target_set.get_html(columns=["name"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        target_set.get_html(columns=["xxx"])
    with pytest.raises(
        ValueError, match=f"Column xxx is not known and cannot be generated."
    ):
        target_set.get_html(columns=["name", "xxx"])
