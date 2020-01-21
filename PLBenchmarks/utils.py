"""
utils.py
Contains utility functions
"""

import numpy as np
import urllib
import json
from simtk import unit



def findPdbUrl(pdb):
    """
    Finds the links to a pdb or a list of pdb codes.

    :param pdb: string
    :return: string compiled string including the urls to the pdb entries
    """

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
    response = urllib.request.urlopen(request)
    page = response.read()
    page = page.decode("utf-8").split()
    res = []
    pdbs = pdb.split()
    for p in page:
        res.append('REP1http://www.rcsb.org/structure/{}REP2{}REP3'.format(p, p))
    for p in pdbs:
        if p not in page:
            raise ValueError(f'PDB {p} not found')
    return ('\n').join(res)



def findDoiUrl(doi):
    """
    Finds the links to a digital object identifier (doi).

    :param doi: string
    :return: string compiled string including the urls to the publication
    """

    url = 'https://api.crossref.org/works/' + str(doi)
    request = urllib.request.Request(url)
    try:
        response = urllib.request.urlopen(request)
        page = response.read().decode("utf-8")
        obj = json.loads(page) 
        if obj['status'] == 'ok':
            obj = obj['message']
        aut = obj['author']
        if len(aut) > 0:
            aut = obj['author'][0]['family']
        else:
            aut = ''
        tit = obj['short-container-title']
        if len(tit) > 0:
            tit = tit[0]
        else:
            tit = ''
        if 'published-print' in obj.keys():
            dat = obj['published-print']['date-parts'][0][0]
        else:
            dat = 'XXXX'
        desc_string='{} et al., {} {}'.format(aut, tit, dat)#, obj['journal-issue']['published-online']['date-parts'][0][0])
        result = f'REP1{obj["URL"]}REP2{desc_string}REP3'
    except urllib.error.HTTPError as e:
        result = doi
    return result


def convertValue(val, originalObs, finalObs, temperature=300.0, outUnit=None):
    """
    Converts an experimental value into another derived quantity with specified unit.
 
    :param val: float, numerical value
    :param originalObs: string, code for the original observable. Can be `dg`, `ki`, `ic50`, `pic50`
    :param finalObs: string, code for the desired derived quantity. Can be `dg`, `ki`, `ic50`, `pic50`
    :param temperature: float, temperature in kelvin
    :param outUnit: type :func:`simtk.unit`, output unit of finalObs, needs to fit to the requested finalObs
    :return: :func:`simtk.unit.quantity` with desired unit
    """

    # define default units
    if outUnit is None:
        if finalObs == 'dg':
            outUnit = unit.kilocalories_per_mole
        elif finalObs == 'ki':
            outUnit = unit.nano * unit.molar
        elif finalObs == 'ic50':
            outUnit = unit.nano * unit.molar
        elif finalObs == 'pic50':
            outUnit = unit.dimensionless
            
    if originalObs == 'dg':
        if finalObs == 'dg':
            return val.in_units_of(outUnit)
        elif finalObs == 'ki':
            return unit.Quantity(np.exp(-val/(unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin))), unit.molar).in_units_of(outUnit)
        elif finalObs == 'ic50':
            return unit.Quantity(np.exp(-val/(unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin))), unit.molar).in_units_of(outUnit)
        elif finalObs == 'pic50':
            return unit.Quantity(val/(unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin))/np.log(10), unit.dimensionless)
        else:
            raise NotImplementedError
    elif originalObs == 'ki':
        if finalObs == 'dg':
            if val.value_in_unit(unit.molar) < 1e-15:
                return unit.Quantity(0.0, outUnit)
            else:
                return (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin) * np.log(val.value_in_unit(unit.molar))).in_units_of(outUnit)
        elif finalObs == 'ki':
            return val.in_units_of(outUnit)
        elif finalObs == 'ic50':
            return val.in_units_of(outUnit)
        elif finalObs == 'pic50':
            if val.value_in_unit(unit.molar) < 1e-15:
                return unit.Quantity(-1e15, outUnit)
            else:
                return unit.Quantity(-np.log(val.value_in_unit(unit.molar))/np.log(10), unit.dimensionless)
        else:
            raise NotImplementedError 
    elif originalObs == 'ic50':
        if finalObs == 'dg':
            if val.value_in_unit(unit.molar) < 1e-15:
                return unit.Quantity(0.0, outUnit)
            else:
                return (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin) * np.log(val.value_in_unit(unit.molar))).in_units_of(outUnit)
        elif finalObs == 'ki':
            return val.in_units_of(outUnit)
        elif finalObs == 'ic50':
            return val.in_units_of(outUnit)
        elif finalObs == 'pic50':
            if val.value_in_unit(unit.molar) < 1e-15:
                return unit.Quantity(-1e15, outUnit)
            else:
                return unit.Quantity(-np.log(val.value_in_unit(unit.molar))/np.log(10), unit.dimensionless)
        else:
            raise NotImplementedError
    elif originalObs == 'pic50':
        if finalObs == 'dg':
            return (-unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(temperature, unit.kelvin) * val.value_in_unit(unit.dimensionless)*np.log(10)).in_units_of(outUnit)
        elif finalObs == 'ki':
            return unit.Quantity(10**(-val.value_in_unit(unit.dimensionless)), unit.molar).in_units_of(outUnit)
        elif finalObs == 'ic50':
            return unit.Quantity(10**(-val.value_in_unit(unit.dimensionless)), unit.molar).in_units_of(outUnit)
        elif finalObs == 'pic50':
            return val
        else:
            raise NotImplementedError     
