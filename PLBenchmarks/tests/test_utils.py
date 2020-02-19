"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import utils
import pytest


def test_findPdbUrl():
    """tests the findPdbUrl function"""
    pdbs = ["6R8L", "6N5O", "6SVL", "6J5R", "6J6D"]
    for pdb in pdbs:
        assert (
            f"REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3"
            == utils.findPdbUrl(pdb)
        )

    strings = [f"REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3" for pdb in pdbs]
    for string in strings:
        assert string in utils.findPdbUrl(" ".join(pdbs)).split()

    with pytest.raises(ValueError, match=f"PDB fakepdb not found"):
        utils.findPdbUrl("fakepdb")

    pdbs[3] = "fakepdb2"
    with pytest.raises(ValueError):
        utils.findPdbUrl(" ".join(pdbs))


def test_findDoiUrl():
    assert (
        "REP1http://dx.doi.org/10.1021/acs.jctc.8b00640REP2Mobley et al., J. Chem. Theory Comput. 2018REP3"
        == utils.findDoiUrl("10.1021/acs.jctc.8b00640")
    )

    assert "fakeDOI" == utils.findDoiUrl("fakeDOI")


def test_convertValue():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'dg'
    ##############################################
    dg = utils.ureg.Quantity(1, utils.ureg("kJ / mole"))
    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(dg, "dg", "dg", outUnit=utils.ureg("kJ / mole")).magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(
            dg, "dg", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.239, eps)
        == utils.convertValue(
            dg, "dg", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.6697, eps)
        == utils.convertValue(dg, "dg", "ki").to(utils.ureg.molar).magnitude
    )
    assert (
        pytest.approx(0.6697, eps)
        == utils.convertValue(dg, "dg", "ic50").to(utils.ureg.molar).magnitude
    )
    assert pytest.approx(0.1741, eps) == utils.convertValue(dg, "dg", "pic50").magnitude
    with pytest.raises(NotImplementedError):
        assert "0.24 kcal/mol" == utils.convertValue(dg, "dg", "fakeObs")

    ##############################################
    # ORIGINAL = 'ki'
    ##############################################
    ki = utils.ureg.Quantity(1, utils.ureg.molar)
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            ki, "ki", "dg", temperature=300, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            ki, "ki", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(ki, "ki", "dg", outUnit=utils.ureg("kcal / mole"))
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(ki, "ki", "ki", outUnit=utils.ureg.molar)
        .to(utils.ureg.molar)
        .magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(ki, "ki", "ic50").to(utils.ureg.molar).magnitude
    )

    assert pytest.approx(0.0, eps) == utils.convertValue(ki, "ki", "pic50").magnitude

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convertValue(ki, "ki", "fakeObs")

    ki = utils.ureg.Quantity(1, utils.ureg("nanomolar"))
    assert (
        pytest.approx(-51.69, eps)
        == utils.convertValue(
            ki, "ki", "dg", temperature=300, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(ki, "ki", "dg").to(utils.ureg("kcal / mole")).magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(ki, "ki", "dg", temperature=300)
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(
            ki, "ki", "dg", temperature=300, outUnit=utils.ureg("kcal / mole")
        )
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-47.04, eps)
        == utils.convertValue(
            ki, "ki", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ki, "ki", "ki").to(utils.ureg("nanomolar")).magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ki, "ki", "ki", outUnit=utils.ureg("nanomolar"))
        .to(utils.ureg("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(ki, "ki", "ic50").to(utils.ureg("nanomolar")).magnitude
    )

    assert pytest.approx(9, eps) == utils.convertValue(ki, "ki", "pic50").magnitude

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convertValue(ki, "ki", "fakeObs")

    ##############################################
    # ORIGINAL = 'ic50'
    ##############################################
    ic50 = utils.ureg.Quantity(1, utils.ureg.molar)
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            ic50, "ic50", "dg", temperature=300, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            ic50, "ic50", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(ic50, "ic50", "dg", outUnit=utils.ureg("kcal / mole"))
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(ic50, "ic50", "ki", outUnit=utils.ureg.molar)
        .to(utils.ureg.molar)
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(ic50, "ic50", "ic50", outUnit=utils.ureg("nanomolar"))
        .to(utils.ureg.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.0, eps) == utils.convertValue(ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convertValue(ic50, "ic50", "fakeObs")

    ic50 = utils.ureg.Quantity(1, utils.ureg("nanomolar"))
    assert (
        pytest.approx(-51.69, eps)
        == utils.convertValue(
            ic50, "ic50", "dg", temperature=300, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(ic50, "ic50", "dg")
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(ic50, "ic50", "dg", temperature=300)
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(
            ic50, "ic50", "dg", temperature=300, outUnit=utils.ureg("kcal / mole")
        )
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-47.04, eps)
        == utils.convertValue(
            ic50, "ic50", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ic50, "ic50", "ki").to(utils.ureg("nanomolar")).magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ic50, "ic50", "ki", outUnit=utils.ureg("nanomolar"))
        .to(utils.ureg("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ic50, "ic50", "ic50")
        .to(utils.ureg("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(ic50, "ic50", "ic50", outUnit=utils.ureg("nanomolar"))
        .to(utils.ureg("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(9.00, eps) == utils.convertValue(ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convertValue(ic50, "ic50", "fakeObs")

    ##############################################
    # ORIGINAL = 'pic50'
    ##############################################
    pic50 = utils.ureg.Quantity(0, "")
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            pic50, "pic50", "dg", temperature=300.0, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(
            pic50, "pic50", "dg", temperature=273.0, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convertValue(pic50, "pic50", "dg", outUnit=utils.ureg("kcal / mole"))
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(pic50, "pic50", "ki", outUnit=utils.ureg.molar)
        .to(utils.ureg.molar)
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convertValue(pic50, "pic50", "ic50", outUnit=utils.ureg.molar)
        .to(utils.ureg.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.0, eps) == utils.convertValue(pic50, "pic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convertValue(pic50, "pic50", "fakeObs")

    pic50 = utils.ureg.Quantity(9, "")
    assert (
        pytest.approx(-51.69, eps)
        == utils.convertValue(
            pic50, "pic50", "dg", temperature=300, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(pic50, "pic50", "dg")
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(pic50, "pic50", "dg", temperature=300)
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convertValue(
            pic50, "pic50", "dg", temperature=300, outUnit=utils.ureg("kcal / mole")
        )
        .to(utils.ureg("kcal / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-47.04, eps)
        == utils.convertValue(
            pic50, "pic50", "dg", temperature=273, outUnit=utils.ureg("kJ / mole")
        )
        .to(utils.ureg("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(pic50, "pic50", "ki")
        .to(utils.ureg("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(pic50, "pic50", "ic50")
        .to(utils.ureg("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convertValue(pic50, "pic50", "ic50", outUnit=utils.ureg("nanomolar"))
        .to(utils.ureg("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(9, eps) == utils.convertValue(pic50, "pic50", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        utils.convertValue(pic50, "pic50", "fakeObs")
