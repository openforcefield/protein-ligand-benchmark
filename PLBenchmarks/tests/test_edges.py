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
    edg = edges.edge(test_dict)
    assert edg.getName() == f"edge_xxx_yyy"
    assert edg.getDict() == {f"edge_xxx_yyy": ["lig_xxx", "lig_yyy"]}

    ligs = ligands.ligandSet("mcl1")
    test_dict = ["50", "60"]
    edg = edges.edge(test_dict)
    assert edg.getName() == f"edge_50_60"
    assert edg.getDict() == {f"edge_50_60": ["lig_50", "lig_60"]}
    pd.testing.assert_series_equal(edg.getDF(), pd.Series({0: "50", 1: "60"}))
    edg.addLigData(ligs)
    with pytest.raises(AssertionError):
        pd.testing.assert_series_equal(edg.getDF(), pd.Series({0: "50", 1: "60"}))
    df = edg.getDF(columns = [0, 1, "Smiles1", "Smiles2", "exp. DeltaG [kcal/mol]", "exp. Error [kcal/mol]"])
    ddg = ligs[f"lig_{edg._data[1]}"]._data[('DerivedMeasurement', 'dg')] \
          - ligs[f"lig_{edg._data[0]}"]._data[('DerivedMeasurement', 'dg')]
    e_ddg = np.sqrt(ligs[f"lig_{edg._data[1]}"]._data[('DerivedMeasurement', 'e_dg')]**2 \
          + ligs[f"lig_{edg._data[0]}"]._data[('DerivedMeasurement', 'e_dg')]**2)
    assert pytest.approx(df["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
    assert pytest.approx(df["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude

def test_edgeSet():
    eps = 0.0001
    ligSet = ligands.ligandSet("mcl1")
    edgSet = edges.edgeSet("mcl1")
    for k, edg in edgSet.items():
        assert f"lig_{edg._data[0]}" in ligSet.keys()
        assert f"lig_{edg._data[1]}" in ligSet.keys()
        ddg = ligSet[f"lig_{edg._data[1]}"]._data[('DerivedMeasurement', 'dg')] \
              - ligSet[f"lig_{edg._data[0]}"]._data[('DerivedMeasurement', 'dg')]
        assert pytest.approx(edg._data["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        e_ddg = np.sqrt(ligSet[f"lig_{edg._data[1]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2 \
                    + ligSet[f"lig_{edg._data[0]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2)
        assert pytest.approx(edg._data["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude

    df = edgSet.getDF()
    for k, edg in df.iterrows():
        assert f"lig_{edg[0]}" in ligSet.keys()
        assert f"lig_{edg[1]}" in ligSet.keys()
        ddg = ligSet[f"lig_{edg[1]}"]._data[('DerivedMeasurement', 'dg')] \
             - ligSet[f"lig_{edg[0]}"]._data[('DerivedMeasurement', 'dg')]
        assert pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude

        e_ddg = np.sqrt(ligSet[f"lig_{edg[1]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2 \
                        + ligSet[f"lig_{edg[0]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2)
        assert pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude

        df2 = edgSet.getDF(columns=[0, 1, "Smiles1", "Smiles2", "exp. DeltaG [kcal/mol]", "exp. Error [kcal/mol]"])
    for k, edg in df2.iterrows():
        assert f"lig_{edg[0]}" in ligSet.keys()
        assert f"lig_{edg[1]}" in ligSet.keys()
        ddg = ligSet[f"lig_{edg[1]}"]._data[('DerivedMeasurement', 'dg')] \
              - ligSet[f"lig_{edg[0]}"]._data[('DerivedMeasurement', 'dg')]
        assert pytest.approx(edg["exp. DeltaG [kcal/mol]"].magnitude, eps) == ddg.magnitude
        e_ddg = np.sqrt(ligSet[f"lig_{edg[1]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2 \
                    + ligSet[f"lig_{edg[0]}"]._data[('DerivedMeasurement', 'e_dg')] ** 2)
        assert pytest.approx(edg["exp. Error [kcal/mol]"].magnitude, 0.5) == e_ddg.magnitude

    html = edgSet.getHTML()
    html = edgSet.getHTML(columns = [0, 1, "Smiles1", "Smiles2", "exp. DeltaG [kcal/mol]", "exp. Error [kcal/mol]"])

    d = edgSet.getDict()
    for edg, ligs in d.items():
        assert edg == f'edge_{ligs[0].replace("lig_", "")}_{ligs[1].replace("lig_", "")}'
        assert ligs[0] in ligSet.keys()
        assert ligs[1] in ligSet.keys()
