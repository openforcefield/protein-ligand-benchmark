"""
utils.py
Contains utility functions
"""

import numpy as np
from scipy import constants
import urllib
import json
from pint import UnitRegistry

import warnings

unit_registry = UnitRegistry()

boltzmann_constant = constants.gas_constant * unit_registry("J / mole / K")


def find_pdb_url(pdb):
    """
    Finds the links to a pdb or a list of pdb codes.

    :param pdb: string
    :return: string compiled string including the urls to the pdb entries
    """

    if pdb is None:
        return ""
    url = "http://www.rcsb.org/pdb/rest/search"
    query_text = f'\
<orgPdbCompositeQuery version="1.0">\
 <queryRefinement>\
   <queryRefinementLevel>0</queryRefinementLevel>\
      <orgPdbQuery>\
        <version>head</version>\
        <queryType>org.pdb.query.simple.StructureIdQuery</queryType>\
        <structureIdList>{pdb}</structureIdList>\
      </orgPdbQuery>\
 </queryRefinement>\
</orgPdbCompositeQuery>\
'
    request = urllib.request.Request(url, data=query_text.encode())
    try:
        response = urllib.request.urlopen(request)
        page = response.read()
        page = page.decode("utf-8").split()
        res = []
        pdb_codes = pdb.split()
        for p in page:
            res.append("REP1http://www.rcsb.org/structure/{}REP2{}REP3".format(p, p))
        for p in pdb_codes:
            if p not in page:
                warnings.warn(f"PDB {p} not found")
    except urllib.error.URLError as e:
        warnings.warn(f"Could not find PDB {pdb}\n{e}")
        res = pdb.split()
    return ("\n").join(res)


def find_doi_url(doi):
    """
    Finds the links to a digital object identifier (doi).

    :param doi: string
    :return: string compiled string including the urls to the publication
    """

    url = "https://api.crossref.org/works/" + str(doi)
    request = urllib.request.Request(url)
    try:
        response = urllib.request.urlopen(request)
        page = response.read().decode("utf-8")
        obj = json.loads(page)
        if obj["status"] == "ok":
            obj = obj["message"]
        aut = obj["author"]
        if len(aut) > 0:
            aut = obj["author"][0]["family"]
        else:
            aut = ""
        tit = obj["short-container-title"]
        if len(tit) > 0:
            tit = tit[0]
        else:
            tit = ""
        if "published-print" in obj.keys():
            dat = obj["published-print"]["date-parts"][0][0]
        else:
            dat = "XXXX"
        desc_string = "{} et al., {} {}".format(
            aut, tit, dat
        )  # , obj['journal-issue']['published-online']['date-parts'][0][0])
        result = f'REP1{obj["URL"]}REP2{desc_string}REP3'
    except urllib.error.URLError as e:
        warnings.warn(f"Could not find DOI: {doi}\n{e}")
        result = doi
    return result


def convert_value(value, original_type, final_type, temperature=300.0, out_unit=None):
    """
    Converts an experimental value into another derived quantity with specified unit.

    :param value: float, numerical value
    :param original_type: string, code for the original observable. Can be `dg`, `ki`, `ic50`, `pic50`
    :param final_type: string, code for the desired derived quantity. Can be `dg`, `ki`, `ic50`, `pic50`
    :param temperature: float, temperature in kelvin
    :param out_unit: unit of type :py:class:`pint`, output unit of final_type, needs to fit to the requested final_type
    :return: :py:class:`pint.Quantity` with desired unit
    """

    # define default units
    if out_unit is None:
        if final_type == "dg":
            out_unit = unit_registry("kilocalories / mole")
        elif final_type == "ki":
            out_unit = unit_registry("nanomolar")
        elif final_type == "ic50":
            out_unit = unit_registry("nanomolar")
        elif final_type == "pic50":
            out_unit = unit_registry("")

    if original_type == "dg":
        if final_type == "dg":
            return value.to(out_unit)
        elif final_type == "ki":
            result = (
                np.exp(
                    -value / (boltzmann_constant * temperature * unit_registry.kelvin)
                )
                * unit_registry.molar
            )
            return result.to(out_unit)
        elif final_type == "ic50":
            result = (
                np.exp(
                    -value / (boltzmann_constant * temperature * unit_registry.kelvin)
                )
                * unit_registry.molar
            )
            return result.to(out_unit)
        elif final_type == "pic50":
            result = (
                value
                / (boltzmann_constant * temperature * unit_registry.kelvin)
                / np.log(10)
            )
            return result.to(out_unit)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ki":
        if final_type == "dg":
            if value < 1e-15 * unit_registry("molar"):
                return 0.0 * out_unit
            else:
                result = (
                    boltzmann_constant
                    * temperature
                    * unit_registry.kelvin
                    * np.log(value / unit_registry.molar)
                )
                return result.to(out_unit).round(2)
        elif final_type == "ki":
            return value.to(out_unit)
        elif final_type == "ic50":
            return value.to(out_unit)
        elif final_type == "pic50":
            if value < 1e-15 * unit_registry("molar"):
                return -1e15 * out_unit
            else:
                result = -np.log(value / unit_registry.molar) / np.log(10)
                return result
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ic50":
        if final_type == "dg":
            if value < 1e-15 * unit_registry("molar"):
                return 0.0 * out_unit
            else:
                result = (
                    boltzmann_constant
                    * temperature
                    * unit_registry.kelvin
                    * np.log(value.to("molar") / unit_registry.molar)
                )
                return result.to(out_unit).round(2)
        elif final_type == "ki":
            return value.to(out_unit)
        elif final_type == "ic50":
            return value.to(out_unit)
        elif final_type == "pic50":
            if value.to("molar") < 1e-15 * unit_registry("molar"):
                return -1e15 * out_unit
            else:
                result = -np.log(value / unit_registry.molar) / np.log(10)
                return result
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "pic50":
        if final_type == "dg":
            result = (
                -boltzmann_constant
                * temperature
                * unit_registry.kelvin
                * value
                * np.log(10)
            )
            return result.to(out_unit).round(2)
        elif final_type == "ki":
            result = 10 ** (-value) * unit_registry("molar")
            return result.to(out_unit)
        elif final_type == "ic50":
            result = 10 ** (-value) * unit_registry("molar")
            return result.to(out_unit)
        elif final_type == "pic50":
            return value.to(out_unit)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )


def convert_error(
    error_value, value, original_type, final_type, temperature=300.0, out_unit=None
):
    """
    Converts an experimental value into another derived quantity with specified unit.

    :param error_value: float, error of val, numerical value
    :param value: float, numerical value
    :param original_type: string, code for the original observable. Can be `dg`, `ki`, `ic50`, `pic50`
    :param final_type: string, code for the desired derived quantity. Can be `dg`, `ki`, `ic50`, `pic50`
    :param temperature: float, temperature in kelvin
    :param out_unit: unit of type :py:class:`pint`, output unit of final_type, needs to fit to the requested final_type
    :return: :py:class:`pint.Quantity` with desired unit
    """

    # define default units
    if out_unit is None:
        if final_type == "dg":
            out_unit = unit_registry("kilocalories / mole")
        elif final_type == "ki":
            out_unit = unit_registry("nanomolar")
        elif final_type == "ic50":
            out_unit = unit_registry("nanomolar")
        elif final_type == "pic50":
            out_unit = unit_registry("")

    if original_type == "dg":
        if final_type == "dg":
            return error_value.to(out_unit)
        elif final_type == "ki":
            # e_ki^2 = (del K/del dG)^2 * e_dG^2
            # e_ki   = 1/RT * exp(-dG/RT) * e_dG
            k_bt = boltzmann_constant * temperature * unit_registry.kelvin
            error = (
                1.0 / k_bt * np.exp(-value / k_bt) * error_value * unit_registry.molar
            )
            return error.to(out_unit)
        elif final_type == "ic50":
            k_bt = boltzmann_constant * temperature * unit_registry.kelvin
            error = (
                1.0 / k_bt * np.exp(-value / k_bt) * error_value * unit_registry.molar
            )
            return error.to(out_unit)
        elif final_type == "pic50":
            # e_pic50^2 = (del pic50/del dG)^2 * e_dG^2
            # e_pic50   = 1/(RT*ln(10)) * e_dG
            k_bt = boltzmann_constant * temperature * unit_registry.kelvin
            error = 1.0 / (k_bt * np.log(10)) * error_value
            return error.to(out_unit)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ki":
        if final_type == "dg":
            if value < 1e-15 * unit_registry.molar:
                return 0.0 * out_unit
            else:
                error = (
                    boltzmann_constant
                    * temperature
                    * unit_registry.kelvin
                    / value
                    * error_value
                )
                return error.to(out_unit).round(2)
        elif final_type == "ki":
            return error_value.to(out_unit)
        elif final_type == "ic50":
            return error_value.to(out_unit)
        elif final_type == "pic50":
            # e_pic50^2 = (del pic50/del Ki)^2 * e_Ki^2
            # e_pic50   = 1/(Ki*ln(10)) * e_Ki
            if (value * np.log(10)) < 1e-15 * unit_registry("molar"):
                return 1e15 * out_unit
            else:
                result = 1 / (value * np.log(10)) * error_value
                return result.to(out_unit).round(2)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ic50":
        if final_type == "dg":
            if value < 1e-15 * unit_registry.molar:
                return 0.0 * out_unit
            else:
                error = (
                    boltzmann_constant
                    * temperature
                    * unit_registry.kelvin
                    / value
                    * error_value
                )
                return error.to(out_unit).round(2)
        elif final_type == "ki":
            return error_value.to(out_unit)
        elif final_type == "ic50":
            return error_value.to(out_unit)
        elif final_type == "pic50":
            # e_pic50^2 = (del pic50/del IC50)^2 * e_IC50^2
            # e_pic50   = 1/(IC50*ln(10)) * e_IC50
            if (value * np.log(10)) < 1e-15 * unit_registry("molar"):
                return 1e15 * out_unit
            else:
                result = 1 / (value * np.log(10)) * error_value
                return result.to(out_unit).round(2)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "pic50":
        if final_type == "dg":
            error = (
                boltzmann_constant
                * temperature
                * unit_registry.kelvin
                * np.log(10)
                * error_value
            )
            return error.to(out_unit).round(2)
        elif final_type == "ki":
            # Ki = 10^(-pIC50)
            # dKi^2 = (del Ki / del pIC50)^2 * dpIC50^2
            # dKi = ln(10) * 10^(-pIC50) * dpIC50
            error = np.log(10) * 10 ** (-value) * error_value * unit_registry("molar")
            return error.to(out_unit).round(2)
        elif final_type == "ic50":
            # IC50 = 10^(-pIC50)
            # dIC50^2 = (del IC50 / del pIC50)^2 * dpIC50^2
            # dIC50 = ln(10) * 10^(-pIC50) * dpIC50
            error = np.log(10) * 10 ** (-value) * error_value * unit_registry("molar")
            return error.to(out_unit).round(2)
        elif final_type == "pic50":
            return error_value.to(out_unit).round(2)
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "\
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
