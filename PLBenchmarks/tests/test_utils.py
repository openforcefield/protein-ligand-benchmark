"""
Unit and regression test for the PLBenchmarks package.
"""

# Import package, test suite, and other packages as needed
from PLBenchmarks import utils
import pytest
import warnings


def test_find_pdb_url():
    """tests the find_pdb_url function"""
    pdbs = ["6R8L", "6N5O", "6SVL", "6J5R", "6J6D"]
    for pdb in pdbs:
        assert (
            f"REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3"
            == utils.find_pdb_url(pdb)
        )

    strings = [f"REP1http://www.rcsb.org/structure/{pdb}REP2{pdb}REP3" for pdb in pdbs]
    for string in strings:
        assert string in utils.find_pdb_url(" ".join(pdbs)).split()

    with pytest.warns(UserWarning):
        utils.find_pdb_url("fakepdb")

    pdbs[3] = "fakepdb2"
    with pytest.warns(UserWarning):
        utils.find_pdb_url(" ".join(pdbs))


def test_find_doi_url():
    assert (
        "REP1http://dx.doi.org/10.1021/acs.jctc.8b00640REP2Mobley et al., J. Chem. Theory Comput. 2018REP3"
        == utils.find_doi_url("10.1021/acs.jctc.8b00640")
    )

    with pytest.warns(UserWarning):
        assert "fakeDOI" == utils.find_doi_url("fakeDOI")

# ToDo @pytest.mark.parametrize("value,final_type,out_unit", [(1.0, "dg", "kJ / mole"])
#def test_convert_value_from_dg(value, final_type, out_unit):
def test_convert_value_from_dg():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'dg'
    ##############################################
    dg = utils.unit_registry.Quantity(1, utils.unit_registry("kJ / mole"))
    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(
            dg, "dg", "dg", out_unit=utils.unit_registry("kJ / mole")
        ).magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(
            dg, "dg", "dg", temperature=273, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.239, eps)
        == utils.convert_value(
            dg, "dg", "dg", temperature=273, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.6697, eps)
        == utils.convert_value(dg, "dg", "ki").to(utils.unit_registry.molar).magnitude
    )
    assert (
        pytest.approx(0.6697, eps)
        == utils.convert_value(dg, "dg", "ic50").to(utils.unit_registry.molar).magnitude
    )
    assert (
        pytest.approx(0.1741, eps) == utils.convert_value(dg, "dg", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        assert "0.24 kcal/mol" == utils.convert_value(dg, "dg", "fakeObs")


def test_convert_value_from_ki():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'ki'
    ##############################################
    ki = utils.unit_registry.Quantity(1, utils.unit_registry.molar)
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ki, "ki", "dg", temperature=300, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ki, "ki", "dg", temperature=273, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ki, "ki", "dg", out_unit=utils.unit_registry("kcal / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(ki, "ki", "ki", out_unit=utils.unit_registry.molar)
        .to(utils.unit_registry.molar)
        .magnitude
    )
    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(ki, "ki", "ic50").to(utils.unit_registry.molar).magnitude
    )

    assert pytest.approx(0.0, eps) == utils.convert_value(ki, "ki", "pic50").magnitude

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convert_value(ki, "ki", "fakeObs")

    ki = utils.unit_registry.Quantity(1, utils.unit_registry("nanomolar"))
    assert (
        pytest.approx(-51.69, eps)
        == utils.convert_value(
            ki, "ki", "dg", temperature=300, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(ki, "ki", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(ki, "ki", "dg", temperature=300)
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(
            ki, "ki", "dg", temperature=300, out_unit=utils.unit_registry("kcal / mole")
        )
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-47.04, eps)
        == utils.convert_value(
            ki, "ki", "dg", temperature=273, out_unit=utils.unit_registry("kJ / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(ki, "ki", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(
            ki, "ki", "ki", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(ki, "ki", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert pytest.approx(9, eps) == utils.convert_value(ki, "ki", "pic50").magnitude

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convert_value(ki, "ki", "fakeObs")

def test_convert_value_from_ic50():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'ic50'
    ##############################################
    ic50 = utils.unit_registry.Quantity(1, utils.unit_registry.molar)
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ic50,
            "ic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            ic50, "ic50", "dg", out_unit=utils.unit_registry("kcal / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(ic50, "ic50", "ki", out_unit=utils.unit_registry.molar)
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(
            ic50, "ic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.0, eps) == utils.convert_value(ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_value(ic50, "ic50", "fakeObs")

    ic50 = utils.unit_registry.Quantity(1, utils.unit_registry("nanomolar"))
    assert (
        pytest.approx(-51.69, eps)
        == utils.convert_value(
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(ic50, "ic50", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(ic50, "ic50", "dg", temperature=300)
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kcal / mole"),
        )
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-47.04, eps)
        == utils.convert_value(
            ic50,
            "ic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(ic50, "ic50", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(
            ic50, "ic50", "ki", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(ic50, "ic50", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(
            ic50, "ic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(9.00, eps) == utils.convert_value(ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_value(ic50, "ic50", "fakeObs")


def test_convert_value_from_ic50():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'pic50'
    ##############################################
    pic50 = utils.unit_registry.Quantity(0, "")
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            pic50,
            "pic50",
            "dg",
            temperature=300.0,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            pic50,
            "pic50",
            "dg",
            temperature=273.0,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(
            pic50, "pic50", "dg", out_unit=utils.unit_registry("kcal / mole")
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(pic50, "pic50", "ki", out_unit=utils.unit_registry.molar)
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(1.0, eps)
        == utils.convert_value(
            pic50, "pic50", "ic50", out_unit=utils.unit_registry.molar
        )
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.0, eps)
        == utils.convert_value(pic50, "pic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_value(pic50, "pic50", "fakeObs")

    pic50 = utils.unit_registry.Quantity(9, "")
    assert (
        pytest.approx(-51.69, eps)
        == utils.convert_value(
            pic50,
            "pic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(pic50, "pic50", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(pic50, "pic50", "dg", temperature=300)
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(-12.35, eps)
        == utils.convert_value(
            pic50,
            "pic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kcal / mole"),
        )
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )

    assert (
        pytest.approx(-47.04, eps)
        == utils.convert_value(
            pic50,
            "pic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(pic50, "pic50", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(pic50, "pic50", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.00, eps)
        == utils.convert_value(
            pic50, "pic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(9, eps) == utils.convert_value(pic50, "pic50", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        utils.convert_value(pic50, "pic50", "fakeObs")


def test_convert_error_from_dg():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'dg'
    ##############################################
    dg = utils.unit_registry.Quantity(1, utils.unit_registry("kJ / mole"))
    edg = utils.unit_registry.Quantity(0.1, utils.unit_registry("kJ / mole"))
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            edg, dg, "dg", "dg", out_unit=utils.unit_registry("kJ / mole")
        ).magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            edg,
            dg,
            "dg",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.0239, eps)
        == utils.convert_error(
            edg,
            dg,
            "dg",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.026849, eps)
        == utils.convert_error(edg, dg, "dg", "ki")
        .to(utils.unit_registry.molar)
        .magnitude
    )
    assert (
        pytest.approx(0.026849, eps)
        == utils.convert_error(edg, dg, "dg", "ic50")
        .to(utils.unit_registry.molar)
        .magnitude
    )
    assert (
        pytest.approx(0.0174, eps)
        == utils.convert_error(edg, dg, "dg", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        assert "0.24 kcal/mol" == utils.convert_error(edg, dg, "dg", "fakeObs")

def test_convert_error_from_ki():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'ki'
    ##############################################
    ki = utils.unit_registry.Quantity(1, utils.unit_registry.molar)
    eki = utils.unit_registry.Quantity(0.1, utils.unit_registry.molar)
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.23, eps)  # 0.226984758
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)  # 0.059616175
        == utils.convert_error(
            eki, ki, "ki", "dg", out_unit=utils.unit_registry("kcal / mole")
        ).magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ki", out_unit=utils.unit_registry.molar)
        .to(utils.unit_registry.molar)
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ic50")
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.04, eps)  # 0.043429448
        == utils.convert_error(eki, ki, "ki", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convert_error(eki, ki, "ki", "fakeObs")

    ki = utils.unit_registry.Quantity(1, utils.unit_registry("nanomolar"))
    eki = utils.unit_registry.Quantity(0.1, utils.unit_registry("nanomolar"))
    assert (
        pytest.approx(0.25, eps)  # 0.02494338
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eki, ki, "ki", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kcal / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.23, eps)
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eki, ki, "ki", "ki", out_unit=utils.unit_registry("nanomolar")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(0.04, eps)
        == utils.convert_error(eki, ki, "ki", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convert_error(eki, ki, "ki", "fakeObs")


def test_convert_error_from_ic50():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'ic50'
    ##############################################
    ic50 = utils.unit_registry.Quantity(1, utils.unit_registry.molar)
    eic50 = utils.unit_registry.Quantity(0.1, utils.unit_registry.molar)
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.23, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        )
        .to(utils.unit_registry("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "dg", out_unit=utils.unit_registry("kcal / mole")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ki", out_unit=utils.unit_registry.molar
        )
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.04, eps)
        == utils.convert_error(eic50, ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_error(eic50, ic50, "ic50", "fakeObs")

    ic50 = utils.unit_registry.Quantity(1, utils.unit_registry("nanomolar"))
    eic50 = utils.unit_registry.Quantity(0.1, utils.unit_registry("nanomolar"))
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eic50, ic50, "ic50", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eic50, ic50, "ic50", "dg", temperature=300)
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kcal / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(0.23, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eic50, ic50, "ic50", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ki", out_unit=utils.unit_registry("nanomolar")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eic50, ic50, "ic50", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        ).magnitude
    )

    assert (
        pytest.approx(0.04, eps)
        == utils.convert_error(eic50, ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_error(eic50, ic50, "ic50", "fakeObs")


def test_convert_error_from_pic50():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'pic50'
    ##############################################
    pic50 = utils.unit_registry.Quantity(0, "")
    epic50 = utils.unit_registry.Quantity(0.5, "")
    assert (
        pytest.approx(2.87, eps)  # 2.871712748
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=300.0,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(2.61, eps)  # 2.613258601
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=273.0,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(
            epic50, pic50, "pic50", "dg", out_unit=utils.unit_registry("kcal / mole")
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ki", out_unit=utils.unit_registry.molar
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ic50", out_unit=utils.unit_registry.molar
        )
        .to(utils.unit_registry.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.5, eps)
        == utils.convert_error(epic50, pic50, "pic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_error(epic50, pic50, "pic50", "fakeObs")

    pic50 = utils.unit_registry.Quantity(9, "")
    epic50 = utils.unit_registry.Quantity(0.5, "")
    assert (
        pytest.approx(2.87, eps)  # 2.871712748
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(epic50, pic50, "pic50", "dg")
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(epic50, pic50, "pic50", "dg", temperature=300)
        .to(utils.unit_registry("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=300,
            out_unit=utils.unit_registry("kcal / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(2.61, eps)  # 2.613258601
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=273,
            out_unit=utils.unit_registry("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(epic50, pic50, "pic50", "ki")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(epic50, pic50, "pic50", "ic50")
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ic50", out_unit=utils.unit_registry("nanomolar")
        )
        .to(utils.unit_registry("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(0.5, eps)
        == utils.convert_error(epic50, pic50, "pic50", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        utils.convert_error(epic50, pic50, "pic50", "fakeObs")
