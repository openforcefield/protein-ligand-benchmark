"""
Unit and regression test for the plbenchmark package.
"""

# Import package, test suite, and other packages as needed
import os
import pytest

import pandas as pd
import yaml
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS
from openff.toolkit.topology import Molecule

import plbenchmark
from plbenchmark import ligands, targets, utils


def test_affinity_data():
    targets.set_data_dir(os.path.join(plbenchmark.__path__[0], "sample_data"))
    file = open(
        os.path.join(targets.get_target_data_path("mcl1_sample") + "ligands.yml")
    )
    data = yaml.full_load(file)
    dfs = []
    for key, d in data.items():
        lig = ligands.Ligand(d)
        lig.derive_observables(
            derived_type="dg", out_unit=utils.unit("kcal / mole")
        )
        lig.add_mol_to_frame()
        lig.get_image()
        dfs.append(lig.get_dataframe(["name", "ROMol", "DerivedMeasurement"]))
    df = pd.DataFrame(dfs)
    assert df.shape[0] == 15

    for n in df.name:
        assert n in [
            "lig_23",
            "lig_26",
            "lig_27",
            "lig_28",
            "lig_29",
            "lig_30",
            "lig_31",
            "lig_32",
            "lig_33",
            "lig_34",
            "lig_35",
            "lig_36",
            "lig_37",
            "lig_38",
            "lig_39",
            # "lig_40",
            # "lig_41",
            # "lig_42",
            # "lig_43",
            # "lig_44",
            # "lig_45",
            # "lig_46",
            # "lig_47",
            # "lig_48",
            # "lig_49",
            # "lig_50",
            # "lig_51",
            # "lig_52",
            # "lig_53",
            # "lig_54",
            # "lig_56",
            # "lig_57",
            # "lig_58",
            # "lig_60",
            # "lig_61",
            # "lig_62",
            # "lig_63",
            # "lig_64",
            # "lig_65",
            # "lig_66",
            # "lig_67",
            # "lig_68",
        ]

    # Check whether the values in the repo are the same and correctly converted by comparing to the values in the JACS paper
    jacs_data = {
        "lig_23": -8.83,
        "lig_26": -8.24,
        "lig_27": -6.12,
        "lig_28": -6.62,
        "lig_29": -6.94,
        "lig_30": -7.85,
        "lig_31": -7.92,
        "lig_32": -6.58,
        "lig_33": -6.88,
        "lig_34": -6.87,
        "lig_35": -8.81,
        "lig_36": -8.18,
        "lig_37": -8.95,
        "lig_38": -7.02,
        "lig_39": -7.03,
        # "lig_40": -7.25,
        # "lig_41": -7.13,
        # "lig_42": -8.9,
        # "lig_43": -7.03,
        # "lig_44": -8.67,
        # "lig_45": -8.95,
        # "lig_46": -7.6,
        # "lig_47": -5.78,
        # "lig_48": -6.66,
        # "lig_49": -8.36,
        # "lig_50": -9.33,
        # "lig_51": -8.45,
        # "lig_52": -9.23,
        # "lig_53": -9.96,
        # "lig_54": -9.78,
        # "lig_56": -9.26,
        # "lig_57": -9.04,
        # "lig_58": -9.41,
        # "lig_60": -8.92,
        # "lig_61": -8.08,
        # "lig_62": -7.96,
        # "lig_63": -9.06,
        # "lig_64": -9.5,
        # "lig_65": -8.41,
        # "lig_66": -8.43,
        # "lig_67": -7.58,
        # "lig_68": -7.69,
    }

    eps = 0.01
    df.index = df.name
    for key, item in jacs_data.items():
        assert (
            pytest.approx(item, eps)
            == df.loc[key, ("DerivedMeasurement", "value")]
            .to(utils.unit("kcal / mole"))
            .magnitude
        )


test_set = []
ligand_set = ligands.LigandSet("mcl1_sample")
for name, lig in ligand_set.items():
    test_set.append(("mcl1_sample", name, lig))


@pytest.mark.parametrize("target, ligand_name, lig", test_set)
def test_ligand_data(target, ligand_name, lig):
    m1 = Chem.MolFromSmiles(lig._data["smiles"][0])
    m1 = Chem.AddHs(m1)
    m2 = Chem.SDMolSupplier(
        os.path.join(
            targets.data_path,
            targets.get_target_dir(target),
            "02_ligands",
            ligand_name,
            "crd",
            f"{ligand_name}.sdf",
        ),
        removeHs=False,
    )[0]
    assert m1.GetNumAtoms() == m2.GetNumAtoms()
    m1.RemoveAllConformers()
    m2.RemoveAllConformers()
    assert pytest.approx(1.0, 1e-9) == DataStructs.FingerprintSimilarity(
        Chem.RDKFingerprint(m1), Chem.RDKFingerprint(m2)
    )
    #            assert Chem.MolToMolBlock(m1) == Chem.MolToMolBlock(m2)
    res = rdFMCS.FindMCS([m1, m2])
    assert res.numAtoms == m1.GetNumAtoms()
    assert res.numBonds == m1.GetNumBonds()

    m3 = lig.get_molecule()
    m2 = Molecule.from_rdkit(m2)
    assert Molecule.are_isomorphic(m2, m3)


def test_ligand_class():
    for target in targets.target_dict.keys():
        ligand_set = ligands.LigandSet(target)
        for name, lig in ligand_set.items():
            assert lig.get_name() == name
            df = lig.get_dataframe()
            assert df["name"][0] == name
            df = lig.get_dataframe(columns=["name"])
            assert df["name"][0] == name
            # ToDo: make proper tests (?)
            lig.find_links()
            lig.get_image()
            lig.get_html()
            lig.get_html(columns=["name", "smiles"])


def test_derive_observables():
    for target in targets.target_dict.keys():
        ligand_set = ligands.LigandSet(target)
        for name, lig in ligand_set.items():
            for i, t in enumerate(["dg", "ki", "ic50", "pic50"]):
                lig.derive_observables(
                    derived_type=t, destination=f"DerivedMeasurement{i}"
                )
                for original_type in lig._observables:
                    if ("measurement", original_type) in list(lig._data.index):
                        assert lig._data[
                            (f"DerivedMeasurement{i}", t)
                        ] == utils.convert_value(
                            lig._data[("measurement", original_type)],
                            original_type=original_type,
                            final_type=t,
                        )

            # Test expected exception when trying to convert to unknown observable
            with pytest.raises(
                NotImplementedError,
                match=f"Conversion to observable xxx not possible. "
                f"Observable must be any of: dg, ki, ic50 or pic50.",
            ):
                lig.derive_observables(
                    derived_type="xxx", destination=f"DerivedMeasurement"
                )

            # Test expected exception when trying to convert from unknown observable
            for original_type in lig._observables:
                if ("measurement", original_type) in list(lig._data.index):
                    lig._data.rename({original_type: "xxx"}, inplace=True, level=1)
                    with pytest.raises(
                        ValueError,
                        match=f"No known measured observable found. "
                        f"Measured observable should be any of: dg, ki, ic50 or pic50.",
                    ):
                        lig.derive_observables(
                            derived_type="pic50", destination=f"DerivedMeasurement{i}"
                        )


def test_ligand_set():
    ligand_set = ligands.LigandSet("mcl1_sample")

    lig_list = ligand_set.get_list()
    for key in lig_list:
        assert key in ligand_set.keys()
        assert isinstance(ligand_set.get_ligand(key), ligands.Ligand)

    with pytest.raises(ValueError, match="Ligand xxx is not part of set."):
        ligand_set.get_ligand("xxx")

    df = ligand_set.get_dataframe()
    for i, row in df.iterrows():
        test_data = row.loc[ligand_set[row.loc["name"][0]]._data.index]
        pd.testing.assert_series_equal(
            ligand_set[row.loc["name"][0]]._data, test_data, check_names=False
        )

    df = ligand_set.get_dataframe(columns=["name", "smiles"])
    for i, row in df.iterrows():
        assert row["name"][0] == ligand_set[row.loc["name"][0]].get_name()
        assert row["smiles"][0] == ligand_set[row.loc["name"][0]]._data["smiles"][0]

    molecules = ligand_set.get_molecules()
    for name, lig in ligand_set.items():
        assert Molecule.are_isomorphic(lig.get_molecule(), molecules[name])

    # ToDo: proper test for get_html()
    ligand_set.get_html()
    ligand_set.get_html(columns=["name", "smiles"])
