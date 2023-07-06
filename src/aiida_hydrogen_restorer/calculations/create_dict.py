# -*- coding: utf-8 -*-
"""Calculation function to create a dictionary"""

import numpy as np

from aiida.engine import calcfunction
from aiida import orm
import collections
#from collections import Counter
import lowdimfinder
import json


@calcfunction
def create_infodict(
    structure_data: orm.StructureData
    )-> dict:
    """Evaluate the dimensionality of the sub components of the structure"""

    s = structure_data.get_ase()
    l = lowdimfinder.LowDimFinder(s, bond_margin=0.2)
    #It is necessary to store the key of the dimensionality as a string. Thus I have rebuild the dictionary (maybe theres a better way?)
    key_dimensionality_noH = str(list(collections.Counter(l.get_group_data()['dimensionality']).keys()))
    value_dimensionality_noH = list(collections.Counter(l.get_group_data()['dimensionality']).values())
    dimensionality_noH = {key_dimensionality_noH : value_dimensionality_noH}
    chemical_formula_noH = dict(collections.Counter(l.get_group_data()['chemical_formula']))

    return {'dimensionality_noH': orm.Dict(dimensionality_noH),
            'chemical_formula_noH' : orm.Dict(chemical_formula_noH)
    }