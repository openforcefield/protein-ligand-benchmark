"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
import os
import pytest

import pandas as pd
import yaml
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS

import PLBenchmarks
from PLBenchmarks import ligands, targets, utils


def testLigand():
    targets.setDataDir(os.path.join(PLBenchmarks.__path__[0], "sample_data"))
    file = open(os.path.join(targets.getTargetDataPath('mcl1') + 'ligands.yml'))
    data = yaml.full_load_all(file)
    dfs = []
    for d in data:
        l = ligands.ligand(d)
        l.deriveObservables(derivedObs="dg", outUnit=utils.ureg("kcal / mole"))
        l.addMolToFrame()
        l.getImg()
        dfs.append(l.getDF(["name", "ROMol", "DerivedMeasurement"]))
    df = pd.DataFrame(dfs)
    assert df.shape[0] == 42

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
            "lig_40",
            "lig_41",
            "lig_42",
            "lig_43",
            "lig_44",
            "lig_45",
            "lig_46",
            "lig_47",
            "lig_48",
            "lig_49",
            "lig_50",
            "lig_51",
            "lig_52",
            "lig_53",
            "lig_54",
            "lig_56",
            "lig_57",
            "lig_58",
            "lig_60",
            "lig_61",
            "lig_62",
            "lig_63",
            "lig_64",
            "lig_65",
            "lig_66",
            "lig_67",
            "lig_68"
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
        "lig_40": -7.25,
        "lig_41": -7.13,
        "lig_42": -8.9,
        "lig_43": -7.03,
        "lig_44": -8.67,
        "lig_45": -8.95,
        "lig_46": -7.6,
        "lig_47": -5.78,
        "lig_48": -6.66,
        "lig_49": -8.36,
        "lig_50": -9.33,
        "lig_51": -8.45,
        "lig_52": -9.23,
        "lig_53": -9.96,
        "lig_54": -9.78,
        "lig_56": -9.26,
        "lig_57": -9.04,
        "lig_58": -9.41,
        "lig_60": -8.92,
        "lig_61": -8.08,
        "lig_62": -7.96,
        "lig_63": -9.06,
        "lig_64": -9.5,
        "lig_65": -8.41,
        "lig_66": -8.43,
        "lig_67": -7.58,
        "lig_68": -7.69,
    }

    eps = 0.01
    for key, item in jacs_data.items():
        assert (
                pytest.approx(item, eps)
                == df[df.name == key][("DerivedMeasurement", "dg")]
                .values[0]
                .to(utils.ureg("kcal / mole"))
                .magnitude
        )


for target in targets.target_list:
    ligSet = ligands.ligandSet(target["name"]).getDF(
        columns=["name", "smiles", "docked"]
    )
    testSet = []
    for index, lig in ligSet.iterrows():
        testSet.append((target["name"], lig["name"][0], target["dir"], lig))


@pytest.mark.parametrize("target, ligName, targetDir, lig", testSet)
def test_ligandData(target, ligName, targetDir, lig):
    m1 = Chem.MolFromSmiles(lig["smiles"][0])
    m2 = Chem.SDMolSupplier(
        os.path.join(
            targets.dataDir,
            targetDir,
            "02_ligands",
            lig["name"][0],
            "crd",
            f'{lig["name"][0]}.sdf',
        )
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


def test_ligand_class():
    for target in targets.target_list:
        ligSet = ligands.ligandSet(target["name"])
        for name, lig in ligSet.items():
            lig.getImg()


def test_ligandSet():
    ligs = ligands.ligandSet("mcl1")
    ligs.getDF()
    ligs.getHTML()
