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


def conv_val(inval, og_type, final_type, temp, out_unit, conv):
    """
    Helper function to convert value units for testing purposes.
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

    comp_val = conv_val(dg, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_val


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

    comp_val = conv_val(ki, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_val


@pytest.mark.parametrize(
        'ic50_units, exp, og_type, final_type, temp, out_units, conv',
        [('molar', 0.0, 'ic50', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ic50', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.0, 'ic50', 'dg', None, 'kcal / mole', 'kJ / mole'),
         ('molar', 1.0, 'ic50', 'ki', None, 'molar', 'molar'),
         ('molar', 1.0, 'ic50', 'ic50', None, 'molar', 'molar'),
         ('molar', 0.0, 'ic50', 'pic50', None, None, None),
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
    ic50 = utils.unit.Quantity(1, ic50_units)

    comp_val = conv_val(ic50, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_val


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

    comp_val = conv_val(pic50, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_val


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


def conv_err(inerr, inval, og_type, final_type, temp, out_unit, conv):
    """
    Helper function to convert error units for testing purposes.
    """
    if temp is not None:
        retval = utils.convert_error(
                inerr, inval, og_type, final_type, temperature=temp,
                out_unit=out_unit)
    else:
        retval = utils.convert_error(
                inerr, inval, og_type, final_type, out_unit=out_unit)

    if conv is not None:
        return retval.to(utils.unit(conv)).magnitude
    else:
        return retval.magnitude


@pytest.mark.parametrize(
        'exp, og_type, final_type, temp, out_units, conv',
        [(0.1, 'dg', 'dg', None, 'kJ / mole', None),
         (0.1, 'dg', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         (0.0239, 'dg', 'dg', 273, 'kJ / mole', 'kcal / mole'),
         (0.026849, 'dg', 'ki', None, None, 'molar'),
         (0.026849, 'dg', 'ic50', None, None, 'molar'),
         (0.0174, 'dg', 'pic50', None, None, None)])
def test_convert_error_from_dg(exp, og_type, final_type, temp, out_units,
                               conv):
    dg = utils.unit.Quantity(1, utils.unit("kJ / mole"))
    edg = utils.unit.Quantity(0.1, utils.unit("kJ / mole"))

    comp_err = conv_err(edg, dg, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_err


@pytest.mark.parametrize('in_value, in_err, in_units, in_type',
        [(1, 0.1, 'kJ / mole', 'dg'),
         (1, 0.1, 'molar', 'ki'),
         (1, 0.1, 'nanomolar', 'ki'),
         (1, 0.1, 'molar', 'ic50'),
         (1, 0.1, 'nanomolar', 'ic50'),
         (0, 0.5, "", 'pic50'),
         (9, 0.5, "", 'pic50')])
def test_convert_error_from_x_notimpl(in_value, in_err, in_units, in_type):
    val = utils.unit.Quantity(in_value, utils.unit(in_units))
    err = utils.unit.Quantity(in_err, utils.unit(in_units))
    with pytest.raises(NotImplementedError):
        utils.convert_error(err, val, in_type, 'fakeObs')


@pytest.mark.parametrize(
        'ki_units, exp, og_type, final_type, temp, out_units, conv',
        [('molar', 0.25, 'ki', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.23, 'ki', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.06, 'ki', 'dg', None, 'kcal / mole', None),
         ('molar', 0.1, 'ki', 'ki', None, 'molar', 'molar'),
         ('molar', 0.1, 'ki', 'ic50', None, None, 'molar'),
         ('molar', 0.04, 'ki', 'pic50', None, None, None),
         ('nanomolar', 0.25, 'ki', 'dg', 300, 'kJ / mole', None),
         ('nanomolar', 0.06, 'ki', 'dg', None, None, 'kcal / mole'),
         ('nanomolar', 0.06, 'ki', 'dg', 300, 'kcal / mole', None),
         ('nanomolar', 0.23, 'ki', 'dg', 273, 'kJ / mole', None),
         ('nanomolar', 0.1, 'ki', 'ki', None, None, 'nanomolar'),
         ('nanomolar', 0.1, 'ki', 'ki', None, 'nanomolar', None),
         ('nanomolar', 0.1, 'ki', 'ic50', None, None, 'nanomolar'),
         ('nanomolar', 0.04, 'ki', 'pic50', None, None, None)])
def test_convert_error_from_ki(ki_units, exp, og_type, final_type, temp,
                               out_units, conv):
    ki = utils.unit.Quantity(1, utils.unit(ki_units))
    eki = utils.unit.Quantity(0.1, utils.unit(ki_units))

    comp_err = conv_err(eki, ki, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_err


@pytest.mark.parametrize(
        'ic50_units, exp, og_type, final_type, temp, out_units, conv',
        [('molar', 0.25, 'ic50', 'dg', 300, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.23, 'ic50', 'dg', 273, 'kJ / mole', 'kJ / mole'),
         ('molar', 0.06, 'ic50', 'dg', None, 'kcal / mole', None),
         ('molar', 0.1, 'ic50', 'ki', None, 'molar', 'molar'),
         ('molar', 0.1, 'ic50', 'ic50', None, 'nanomolar', 'molar'),
         ('molar', 0.04, 'ic50', 'pic50', None, None, None),
         ('nanomolar', 0.25, 'ic50', 'dg', 300, 'kJ /mole', None),
         ('nanomolar', 0.06, 'ic50', 'dg', None, 'kcal / mole', None),
         ('nanomolar', 0.06, 'ic50', 'dg', 300, None, 'kcal / mole'),
         ('nanomolar', 0.06, 'ic50', 'dg', 300, 'kcal / mole', None),
         ('nanomolar', 0.23, 'ic50', 'dg', 273, 'kJ / mole', None),
         ('nanomolar', 0.1, 'ic50', 'ki', None, None, 'nanomolar'),
         ('nanomolar', 0.1, 'ic50', 'ki', None, 'nanomolar', None),
         ('nanomolar', 0.1, 'ic50', 'ic50', None, None, 'nanomolar'),
         ('nanomolar', 0.1, 'ic50', 'ic50', None, 'nanomolar', None),
         ('nanomolar', 0.04, 'ic50', 'pic50', None, None, None)])
def test_convert_error_from_ic50(ic50_units, exp, og_type, final_type, temp,
                                 out_units, conv):

    ic50 = utils.unit.Quantity(1, utils.unit(ic50_units))
    eic50 = utils.unit.Quantity(0.1, utils.unit(ic50_units))

    comp_err = conv_err(
        eic50, ic50, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_err


@pytest.mark.parametrize(
    'pic50_val, epic50_val, exp, og_type, final_type, temp, out_units, conv',
    [(0, 0.5, 2.87, 'pic50', 'dg', 300, 'kJ / mole', None),
     (0, 0.5, 2.61, 'pic50', 'dg', 273, 'kJ / mole', None),
     (0, 0.5, 0.69, 'pic50', 'dg', None, 'kcal / mole', None),
     (0, 0.5, 1.15, 'pic50', 'ki', None, 'molar', None),
     (0, 0.5, 1.15, 'pic50', 'ic50', None, 'molar', 'molar'),
     (0, 0.5, 0.5, 'pic50', 'pic50', None, None, None),
     (9, 0.5, 2.87, 'pic50', 'dg', 300, 'kJ /mole', None),
     (9, 0.5, 0.69, 'pic50', 'dg', None, 'kcal / mole', None),
     (9, 0.5, 0.69, 'pic50', 'dg', 300, None, 'kcal / mole'),
     (9, 0.5, 0.69, 'pic50', 'dg', 300, 'kcal / mole', None),
     (9, 0.5, 2.61, 'pic50', 'dg', 273, 'kJ / mole', None),
     (9, 0.5, 1.15, 'pic50', 'ki', None, None, 'nanomolar'),
     (9, 0.5, 1.15, 'pic50', 'ic50', None, None, 'nanomolar'),
     (9, 0.5, 1.15, 'pic50', 'ic50', None, 'nanomolar', 'nanomolar'),
     (9, 0.5, 0.5, 'pic50', 'pic50', None, None, None)])
def test_convert_error_from_pic50(pic50_val, epic50_val, exp, og_type,
                                  final_type, temp, out_units, conv):

    pic50 = utils.unit.Quantity(pic50_val, "")
    epic50 = utils.unit.Quantity(epic50_val, "")

    comp_err = conv_err(
        epic50, pic50, og_type, final_type, temp, out_units, conv)
    assert pytest.approx(exp, 0.001) == comp_err
