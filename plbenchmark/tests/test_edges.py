"""
Unit and regression test for the plbenchmark package.
"""

# Import package, test suite, and other packages as needed
import pytest
import numpy as np
import pandas as pd
from plbenchmark import edges, ligands


def test_edge():
    eps = 0.0001
    test_dict = {"name": "edge_xxx_yyy", "ligand_a": "xxx", "ligand_b": "yyy"}
    edg = edges.Edge(test_dict)
    assert edg.get_name() == "edge_xxx_yyy"
    assert edg.get_dict() == {
        "name": "edge_xxx_yyy",
        "ligand_a": "xxx",
        "ligand_b": "yyy",
    }

    ligand_set = ligands.LigandSet("mcl1_sample")
    test_dict = {"name": "edge_30_27", "ligand_a": "lig_30", "ligand_b": "lig_27"}
    edg = edges.Edge(test_dict)
    assert edg.get_name() == "edge_30_27"
    assert edg.get_dict() == {
        "name": "edge_30_27",
        "ligand_a": "lig_30",
        "ligand_b": "lig_27",
    }
    pd.testing.assert_series_equal(
        edg.get_dataframe(),
        pd.Series({"name": "edge_30_27", "ligand_a": "lig_30", "ligand_b": "lig_27"}),
    )
    edg.add_ligand_data(ligand_set)
    with pytest.raises(AssertionError):
        pd.testing.assert_series_equal(
            edg.get_dataframe(),
            pd.Series(
                {"name": "edge_30_27", "ligand_a": "lig_30", "ligand_b": "lig_27"}
            ),
        )
    df = edg.get_dataframe(
        columns=[
            "ligand_a",
            "ligand_b",
            "Smiles1",
            "Smiles2",
            "exp. DeltaG [kcal/mol]",
            "exp. Error [kcal/mol]",
        ]
    )
    ddg = (
        ligand_set[f"{edg._data['ligand_b']}"]._data[("DerivedMeasurement", "value")]
        - ligand_set[f"{edg._data['ligand_a']}"]._data[("DerivedMeasurement", "value")]
    )
    e_ddg = np.sqrt(
        ligand_set[f"{edg._data['ligand_b']}"]._data[("DerivedMeasurement", "error")]
        ** 2
        + ligand_set[f"{edg._data['ligand_a']}"]._data[("DerivedMeasurement", "error")]
        ** 2
    )
    assert pytest.approx(df["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
    assert pytest.approx(df["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude


def test_edge_set():
    eps = 0.01
    lig_set = ligands.LigandSet("mcl1_sample")
    edg_set = edges.EdgeSet("mcl1_sample")
    for k, edg in edg_set.items():
        assert f"{edg._data['ligand_a']}" in lig_set.keys()
        assert f"{edg._data['ligand_b']}" in lig_set.keys()
        ddg = (
            lig_set[f"{edg._data['ligand_b']}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"{edg._data['ligand_a']}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg._data["exp. DeltaG [kcal/mol]"].magnitude, eps)
            == ddg.magnitude
        )
        e_ddg = np.sqrt(
            lig_set[f"{edg._data['ligand_a']}"]._data[("DerivedMeasurement", "error")]
            ** 2
            + lig_set[f"{edg._data['ligand_b']}"]._data[("DerivedMeasurement", "error")]
            ** 2
        )
        assert (
            pytest.approx(edg._data["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

    df = edg_set.get_dataframe()
    for k, edg in df.iterrows():
        assert f"{edg['ligand_a']}" in lig_set.keys()
        assert f"{edg['ligand_b']}" in lig_set.keys()
        ddg = (
            lig_set[f"{edg['ligand_b']}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"{edg['ligand_a']}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        )

        e_ddg = np.sqrt(
            lig_set[f"{edg['ligand_b']}"]._data[("DerivedMeasurement", "error")] ** 2
            + lig_set[f"{edg['ligand_a']}"]._data[("DerivedMeasurement", "error")] ** 2
        )
        assert (
            pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

        df2 = edg_set.get_dataframe(
            columns=[
                "ligand_a",
                "ligand_b",
                "Smiles1",
                "Smiles2",
                "exp. DeltaG [kcal/mol]",
                "exp. Error [kcal/mol]",
            ]
        )
    for k, edg in df2.iterrows():
        assert f"{edg['ligand_a']}" in lig_set.keys()
        assert f"{edg['ligand_b']}" in lig_set.keys()
        ddg = (
            lig_set[f"{edg['ligand_b']}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"{edg['ligand_a']}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        )

        e_ddg = np.sqrt(
            lig_set[f"{edg['ligand_b']}"]._data[("DerivedMeasurement", "error")] ** 2
            + lig_set[f"{edg['ligand_a']}"]._data[("DerivedMeasurement", "error")] ** 2
        )
        assert (
            pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

    html = edg_set.get_html()
    html = edg_set.get_html(
        columns=[
            "ligand_a",
            "ligand_b",
            "Smiles1",
            "Smiles2",
            "exp. DeltaG [kcal/mol]",
            "exp. Error [kcal/mol]",
        ]
    )

    d = edg_set.get_dict()
    for edg, ligs in d.items():
        assert (
            edg
            == f'edge_{ligs["ligand_a"].replace("lig_", "")}_{ligs["ligand_b"].replace("lig_", "")}'
        )
        assert ligs["ligand_a"] in lig_set.keys()
        assert ligs["ligand_b"] in lig_set.keys()
