"""
util.py
Protein-Ligand Benchmark Dataset for testing Parameters and Methods of Free Energy Calculations.
Utility functions
"""

import numpy as np
import urllib
import json
from simtk import unit

import pandas as pd

def findPdbUrl(pdb):
    """
    Finds the links to a pdb or a list of pdb codes.
    Parameters
    ----------
    pdb: str
    Returns
    -------
    links : str
        Compiled string including the urls to the pdb entries

    """

    url = "http://www.rcsb.org/pdb/rest/search"
    query_text = f"""
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
   <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.StructureIdQuery</queryType>
        <structureIdList>{pdb}</structureIdList>
      </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
"""
    request = urllib.request.Request(url, data=query_text.encode())
    try:
        response = urllib.request.urlopen(request)
        page = response.read()
        #plTable['PDB'].at[i]=page.decode("utf-8")
        page = page.decode("utf-8").split()
        res = []
        for p in page:
             res.append('REP1http://www.rcsb.org/structure/{}REP2{}REP3'.format(p, p))
    except urllib.error.HTTPError as e:
        print(f"PDB error...")
    return ('\n').join(res)



def findDoiUrl(doi):
    """
    Finds the links to a digital object identifier (doi).
    Parameters
    ----------
    doi: str
    Returns
    -------
    link : str
        Compiled string including the urls to the publication

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
        request = urllib.request.Request(str(d))
        try:
            response = urllib.request.urlopen(request)
            result = f'REP1{d}REP2{d}REP3'
        except urllib.error.HTTPError as e:
            print(f"HTTPError...")
            result = d
    return result


def convertValue(val, originalObs, finalObs, energyUnit=unit.kilocalories_per_mole):
    """
    Converts an experimental value into another derived quantity with specified unit.
    Parameters
    ----------
    val: float, numerical value
    originalObs: str, code for the original observable. Can be 'dg', 'ki', 'ic50', 'pic50'
    finalObs: str, code for the desired derived quantity. Can be 'dg', 'ki', 'ic50', 'pic50'
    energyUnit: type simtk.unit
    Returns
    -------
    quantity : simtk.quantity with desired unit

    """
    if originalObs == 'dg':
        if finalObs == 'dg':
            return val.in_units_of(energyUnit).format('%.2f')
        elif finalObs == 'ki':
            raise NotImplementedError
        elif finalObs == 'ic50':
            raise NotImplementedError
        elif finalObs == 'pic50':
            raise NotImplementedError
        else:
            raise NotImplementedError
    elif originalObs == 'ki':
        if finalObs == 'dg':
            return (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(300.0, unit.kelvin) * np.log(val.value_in_unit(unit.molar))).in_units_of(energyUnit).format('%.2f')
        elif finalObs == 'ki':
            return val.in_units_of(unit.nano * unit.molar).format('%.2f')
        elif finalObs == 'ic50':
            raise NotImplementedError
        elif finalObs == 'pic50':
            raise NotImplementedError   
        else:
            raise NotImplementedError 
    elif originalObs == 'ic50':
        if finalObs == 'dg':
            return (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(300.0, unit.kelvin) * np.log(val.value_in_unit(unit.molar))).in_units_of(energyUnit).format('%.2f')
        elif finalObs == 'ki':
            raise NotImplementedError    
        elif finalObs == 'ic50':
            return val.in_units_of(unit.nano * unit.molar).format('%.2f')
        elif finalObs == 'pic50':
            raise NotImplementedError    
        else:
            raise NotImplementedError
    elif originalObs == 'pic50':
        if finalObs == 'dg':
            return (-unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * unit.Quantity(300.0, unit.kelvin) * val.value_in_unit(unit.dimensionless)*np.log(10)).in_units_of(energyUnit).format('%.2f')
        elif finalObs == 'ki':
            raise NotImplementedError
        elif finalObs == 'ic50':
            raise NotImplementedError
        elif finalObs == 'pic50':
            return val.in_units_of(unit.nano * unit.molar).format('%.2f')     
        else:
            raise NotImplementedError     
