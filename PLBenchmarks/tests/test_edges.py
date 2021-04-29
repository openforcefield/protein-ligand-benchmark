"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
import pytest
import numpy as np
import pandas as pd
from PLBenchmarks import edges, ligands


def test_edge():
    eps = 0.0001
    test_dict = ["xxx", "yyy"]
    edg = edges.Edge(test_dict)
    assert edg.get_name() == f"edge_xxx_yyy"
    assert edg.get_dict() == {f"edge_xxx_yyy": ["lig_xxx", "lig_yyy"]}

    ligand_set = ligands.LigandSet("mcl1_sample")
    test_dict = ["30", "27"]
    edg = edges.Edge(test_dict)
    assert edg.get_name() == f"edge_30_27"
    assert edg.get_dict() == {f"edge_30_27": ["lig_30", "lig_27"]}
    pd.testing.assert_series_equal(edg.get_dataframe(), pd.Series({0: "30", 1: "27"}))
    edg.add_ligand_data(ligand_set)
    with pytest.raises(AssertionError):
        pd.testing.assert_series_equal(
            edg.get_dataframe(), pd.Series({0: "30", 1: "27"})
        )
    df = edg.get_dataframe(
        columns=[
            0,
            1,
            "Smiles1",
            "Smiles2",
            "exp. DeltaG [kcal/mol]",
            "exp. Error [kcal/mol]",
        ]
    )
    ddg = (
        ligand_set[f"lig_{edg._data[1]}"]._data[("DerivedMeasurement", "value")]
        - ligand_set[f"lig_{edg._data[0]}"]._data[("DerivedMeasurement", "value")]
    )
    e_ddg = np.sqrt(
        ligand_set[f"lig_{edg._data[1]}"]._data[("DerivedMeasurement", "error")] ** 2
        + ligand_set[f"lig_{edg._data[0]}"]._data[("DerivedMeasurement", "error")] ** 2
    )
    assert pytest.approx(df["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
    assert pytest.approx(df["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude


def test_edge_set():
    eps = 0.01
    lig_set = ligands.LigandSet("mcl1_sample")
    edg_set = edges.EdgeSet("mcl1_sample")
    for k, edg in edg_set.items():
        assert f"lig_{edg._data[0]}" in lig_set.keys()
        assert f"lig_{edg._data[1]}" in lig_set.keys()
        ddg = (
            lig_set[f"lig_{edg._data[1]}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"lig_{edg._data[0]}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg._data["exp. DeltaG [kcal/mol]"].magnitude, eps)
            == ddg.magnitude
        )
        e_ddg = np.sqrt(
            lig_set[f"lig_{edg._data[1]}"]._data[("DerivedMeasurement", "error")] ** 2
            + lig_set[f"lig_{edg._data[0]}"]._data[("DerivedMeasurement", "error")] ** 2
        )
        assert (
            pytest.approx(edg._data["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

    df = edg_set.get_dataframe()
    for k, edg in df.iterrows():
        assert f"lig_{edg[0]}" in lig_set.keys()
        assert f"lig_{edg[1]}" in lig_set.keys()
        ddg = (
            lig_set[f"lig_{edg[1]}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"lig_{edg[0]}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        )

        e_ddg = np.sqrt(
            lig_set[f"lig_{edg[1]}"]._data[("DerivedMeasurement", "error")] ** 2
            + lig_set[f"lig_{edg[0]}"]._data[("DerivedMeasurement", "error")] ** 2
        )
        assert (
            pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

        df2 = edg_set.get_dataframe(
            columns=[
                0,
                1,
                "Smiles1",
                "Smiles2",
                "exp. DeltaG [kcal/mol]",
                "exp. Error [kcal/mol]",
            ]
        )
    for k, edg in df2.iterrows():
        assert f"lig_{edg[0]}" in lig_set.keys()
        assert f"lig_{edg[1]}" in lig_set.keys()
        ddg = (
            lig_set[f"lig_{edg[1]}"]._data[("DerivedMeasurement", "value")]
            - lig_set[f"lig_{edg[0]}"]._data[("DerivedMeasurement", "value")]
        )
        assert (
            pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        )
        e_ddg = np.sqrt(
            lig_set[f"lig_{edg[1]}"]._data[("DerivedMeasurement", "error")] ** 2
            + lig_set[f"lig_{edg[0]}"]._data[("DerivedMeasurement", "error")] ** 2
        )
        assert (
            pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5)
            == e_ddg.magnitude
        )

    html = edg_set.get_html()
    html = edg_set.get_html(
        columns=[
            0,
            1,
            "Smiles1",
            "Smiles2",
            "exp. DeltaG [kcal/mol]",
            "exp. Error [kcal/mol]",
        ]
    )

    d = edg_set.get_dict()
    for edg, ligs in d.items():
        assert (
            edg == f'edge_{ligs[0].replace("lig_", "")}_{ligs[1].replace("lig_", "")}'
        )
        assert ligs[0] in lig_set.keys()
        assert ligs[1] in lig_set.keys()
