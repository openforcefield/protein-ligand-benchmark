"""
Unit and regression test for the plbenchmark package.
"""

# Import package, test suite, and other packages as needed
from plbenchmark import utils
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
        assert string in utils.find_pdb_url(pdbs)

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


def conv_unit(inval, og_type, final_type, temp, out_unit, conv):
    """
    Helper function to convert unit values for testing purposes.
    """
    if temp is not None:
        retval = utils.convert_value(
                inval, og_type, final_type, temperature=temp,
                out_unit=out_unit)
    else:
        retval = utils.convert_value(
                inval, og_type, final_type, out_unit=out_unit)

    if conv is not None:
        return retval.to(utils.unit(conv)).magnitude
    else:
        return retval.magnitude


@pytest.mark.parametrize("exp, og_type, final_type, temp, out_units, conv", [
    (1.0, "dg", "dg", None, "kJ / mole", None),
    (1.0, "dg", "dg", 273.0, "kJ / mole", "kJ / mole"),
    (0.239, "dg", "dg", 273.0, "kJ / mole", "kcal / mole"),
    (0.6697, "dg", "ki", None, None, "molar"),
    (0.6697, "dg", "ic50", None, None, "molar"),
    (0.1741, "dg", "pic50", None, None, None)])
def test_convert_value_from_dg(exp, og_type, final_type, temp, out_units,
                               conv):
    dg = utils.unit.Quantity(1, utils.unit("kJ / mole"))

    conv_val = conv_unit(dg, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == conv_val


@pytest.mark.parametrize(
        'ki_units, exp, og_type, final_type, temp, out_units, conv',
        [('molar', 0.0, 'ki', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ki', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ki', 'dg', None, 'kcal / mole', 'kJ / mole'),
         ('molar', 1.0, 'ki', 'ki', None, 'molar', 'molar'),
         ('molar', 1.0, 'ki', 'ic50', None, None, 'molar'),
         ('molar', 0.0, 'ki', 'pic50', None, None, None),
         ('nanomolar', -51.69, 'ki', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('nanomolar', -12.35, 'ki', 'dg', None, None, 'kcal / mole'),
         ('nanomolar', -12.35, 'ki', 'dg', 300, None, 'kcal / mole'),
         ('nanomolar', -12.35, 'ki', 'dg', 300, 'kcal / mole', 'kcal / mole'),
         ('nanomolar', -47.04, 'ki', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('nanomolar', 1.0, 'ki', 'ki', None, None, 'nanomolar'),
         ('nanomolar', 1.0, 'ki', 'ki', None, 'nanomolar', 'nanomolar'),
         ('nanomolar', 1.0, 'ki', 'ic50', None, None, 'nanomolar'),
         ('nanomolar', 9.0, 'ki', 'pic50', None, None, None)])
def test_convert_value_from_ki(ki_units, exp, og_type, final_type, temp,
                               out_units, conv):
    ki = utils.unit.Quantity(1, ki_units)

    conv_val = conv_unit(ki, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == conv_val


@pytest.mark.parametrize(
        'ic50_units, exp, og_type, final_type, temp, out_units, conv',
        [('molar', 0.0, 'ic50', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ic50', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ic50', 'dg', None, 'kcal / mole', 'kJ / mole'),
         ('molar', 1.0, 'ic50', 'ki', None, 'molar', 'molar'),
         ('molar', 1.0, 'ic50', 'ic50', None, 'molar', 'molar'),
         ('molar', 1.0, 'ic50', 'pic50', None, None, None),
         ('nanomolar', -51.69, 'ic50', 'dg', 300, 'kJ / mole', 'kJ/mole'),
         ('nanomolar', -12.35, 'ic50', 'dg', None, None, 'kcal / mole'),
         ('nanomolar', -12.35, 'ic50', 'dg', 300, None, 'kcal / mole'),
         ('nanomolar', -12.35, 'ic50', 'dg', 300, 'kcal / mole', 'kcal / mole'),
         ('nanomolar', -47.04, 'ic50', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('nanomolar', 1.0, 'ic50', 'ki', None, None, 'nanomolar'),
         ('nanomolar', 1.0, 'ic50', 'ki', None, 'nanomolar', 'nanomolar'),
         ('nanomolar', 1.0, 'ic50', 'ic50', None, None, 'nanomolar'),
         ('nanomolar', 1.0, 'ic50', 'ic50', None, 'nanomolar', 'nanomolar'),
         ('nanomolar', 9.0, 'ic50', 'pic50', None, None, None)])
def test_convert_value_from_ic50(ic50_units, exp, og_type, final_type, temp,
                                 out_units, conv):
    ic50 = utils.unit.Quantity(1, utils.unit.molar)

    conv_val = conv_unit(ic50, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == conv_val


@pytest.mark.parametrize(
        'pic50_val, exp, og_type, final_type, temp, out_units, conv',
        [(0, 0.0, 'pic50', 'dg', 300.0, 'kJ / mole', 'kJ / mole'),
         (0, 0.0, 'pic50', 'dg', 273.0, 'kJ / mole', 'kJ / mole'),
         (0, 0.0, 'pic50', 'dg', None, 'kcal / mole', 'kJ / mole'),
         (0, 1.0, 'pic50', 'ki', None, 'molar', 'molar'),
         (0, 1.0, 'pic50', 'ic50', None, 'molar', 'molar'),
         (0, 0.0, 'pic50', 'pic50', None, None, None),
         (9, -51.69, 'pic50', 'dg', 300, 'kJ / mole', 'kJ /mole'),
         (9, -12.35, 'pic50', 'dg', None, None, 'kcal /mole'),
         (9, -12.35, 'pic50', 'dg', 300, None, 'kcal /mole'),
         (9, -12.35, 'pic50', 'dg', 300, 'kcal / mole', 'kcal /mole'),
         (9, -47.04, 'pic50', 'dg', 273, 'kJ / mole', 'kJ /mole'),
         (9, 1.0, 'pic50', 'ki', None, None, 'nanomolar'),
         (9, 1.0, 'pic50', 'ic50', None, None, 'nanomolar'),
         (9, 1.0, 'pic50', 'ic50', None, 'nanomolar', 'nanomolar'),
         (9, 9.0, 'pic50', 'pic50', None, None, None)])
def test_convert_value_from_pic50(pic50_val, exp, og_type, final_type, temp,
                                  out_units, conv):
    pic50 = utils.unit.Quantity(pic50_val, "")

    conv_val = conv_unit(pic50, og_type, final_type, temp, out_units, conv)


@pytest.mark.parametrize('value, in_units, in_type',
        [(1, "kJ / mole", "dg"),
         (1, "molar", "ki"),
         (1, "nanomolar", "ki"),
         (1, "molar", "ic50"),
         (1, "nanomolar", "ic50"),
         (0, "", "pic50"),
         (9, "", "pic50")])
def test_convert_value_from_x_not_implemented_error(value, in_units, in_type):
    var = utils.unit.Quantity(value, utils.unit(in_units))

    with pytest.raises(NotImplementedError):
        utils.convert_value(var, in_type, "fakeObs")


def test_convert_error_from_dg():
    eps = 0.001
    ##############################################
    # ORIGINAL = 'dg'
    ##############################################
    dg = utils.unit.Quantity(1, utils.unit("kJ / mole"))
    edg = utils.unit.Quantity(0.1, utils.unit("kJ / mole"))
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            edg, dg, "dg", "dg", out_unit=utils.unit("kJ / mole")
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
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kJ / mole"))
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
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.026849, eps)
        == utils.convert_error(edg, dg, "dg", "ki")
        .to(utils.unit.molar)
        .magnitude
    )
    assert (
        pytest.approx(0.026849, eps)
        == utils.convert_error(edg, dg, "dg", "ic50")
        .to(utils.unit.molar)
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
    ki = utils.unit.Quantity(1, utils.unit.molar)
    eki = utils.unit.Quantity(0.1, utils.unit.molar)
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=300,
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kJ / mole"))
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
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)  # 0.059616175
        == utils.convert_error(
            eki, ki, "ki", "dg", out_unit=utils.unit("kcal / mole")
        ).magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ki", out_unit=utils.unit.molar)
        .to(utils.unit.molar)
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ic50")
        .to(utils.unit.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.04, eps)  # 0.043429448
        == utils.convert_error(eki, ki, "ki", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        assert "xxx" == utils.convert_error(eki, ki, "ki", "fakeObs")

    ki = utils.unit.Quantity(1, utils.unit("nanomolar"))
    eki = utils.unit.Quantity(0.1, utils.unit("nanomolar"))
    assert (
        pytest.approx(0.25, eps)  # 0.02494338
        == utils.convert_error(
            eki,
            ki,
            "ki",
            "dg",
            temperature=300,
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eki, ki, "ki", "dg")
        .to(utils.unit("kcal / mole"))
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
            out_unit=utils.unit("kcal / mole"),
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
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ki")
        .to(utils.unit("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eki, ki, "ki", "ki", out_unit=utils.unit("nanomolar")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eki, ki, "ki", "ic50")
        .to(utils.unit("nanomolar"))
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
    ic50 = utils.unit.Quantity(1, utils.unit.molar)
    eic50 = utils.unit.Quantity(0.1, utils.unit.molar)
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kJ / mole"))
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
            out_unit=utils.unit("kJ / mole"),
        )
        .to(utils.unit("kJ / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "dg", out_unit=utils.unit("kcal / mole")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ki", out_unit=utils.unit.molar
        )
        .to(utils.unit.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ic50", out_unit=utils.unit("nanomolar")
        )
        .to(utils.unit.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.04, eps)
        == utils.convert_error(eic50, ic50, "ic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_error(eic50, ic50, "ic50", "fakeObs")

    ic50 = utils.unit.Quantity(1, utils.unit("nanomolar"))
    eic50 = utils.unit.Quantity(0.1, utils.unit("nanomolar"))
    assert (
        pytest.approx(0.25, eps)
        == utils.convert_error(
            eic50,
            ic50,
            "ic50",
            "dg",
            temperature=300,
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eic50, ic50, "ic50", "dg")
        .to(utils.unit("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.06, eps)
        == utils.convert_error(eic50, ic50, "ic50", "dg", temperature=300)
        .to(utils.unit("kcal / mole"))
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
            out_unit=utils.unit("kcal / mole"),
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
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eic50, ic50, "ic50", "ki")
        .to(utils.unit("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ki", out_unit=utils.unit("nanomolar")
        ).magnitude
    )

    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(eic50, ic50, "ic50", "ic50")
        .to(utils.unit("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(0.1, eps)
        == utils.convert_error(
            eic50, ic50, "ic50", "ic50", out_unit=utils.unit("nanomolar")
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
    pic50 = utils.unit.Quantity(0, "")
    epic50 = utils.unit.Quantity(0.5, "")
    assert (
        pytest.approx(2.87, eps)  # 2.871712748
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=300.0,
            out_unit=utils.unit("kJ / mole"),
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
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )
    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(
            epic50, pic50, "pic50", "dg", out_unit=utils.unit("kcal / mole")
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ki", out_unit=utils.unit.molar
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ic50", out_unit=utils.unit.molar
        )
        .to(utils.unit.molar)
        .magnitude
    )

    assert (
        pytest.approx(0.5, eps)
        == utils.convert_error(epic50, pic50, "pic50", "pic50").magnitude
    )

    with pytest.raises(NotImplementedError):
        utils.convert_error(epic50, pic50, "pic50", "fakeObs")

    pic50 = utils.unit.Quantity(9, "")
    epic50 = utils.unit.Quantity(0.5, "")
    assert (
        pytest.approx(2.87, eps)  # 2.871712748
        == utils.convert_error(
            epic50,
            pic50,
            "pic50",
            "dg",
            temperature=300,
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(epic50, pic50, "pic50", "dg")
        .to(utils.unit("kcal / mole"))
        .magnitude
    )
    assert (
        pytest.approx(0.69, eps)  # 0.686356035
        == utils.convert_error(epic50, pic50, "pic50", "dg", temperature=300)
        .to(utils.unit("kcal / mole"))
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
            out_unit=utils.unit("kcal / mole"),
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
            out_unit=utils.unit("kJ / mole"),
        ).magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(epic50, pic50, "pic50", "ki")
        .to(utils.unit("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(epic50, pic50, "pic50", "ic50")
        .to(utils.unit("nanomolar"))
        .magnitude
    )
    assert (
        pytest.approx(1.15, eps)  # 1.151292546
        == utils.convert_error(
            epic50, pic50, "pic50", "ic50", out_unit=utils.unit("nanomolar")
        )
        .to(utils.unit("nanomolar"))
        .magnitude
    )

    assert (
        pytest.approx(0.5, eps)
        == utils.convert_error(epic50, pic50, "pic50", "pic50").magnitude
    )
    with pytest.raises(NotImplementedError):
        utils.convert_error(epic50, pic50, "pic50", "fakeObs")
