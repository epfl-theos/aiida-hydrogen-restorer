# -*- coding: utf-8 -*-
"""Calculation function to add hydrogen to a structure based on the potential."""

from skimage.feature.peak import peak_local_max
import numpy as np

from aiida.engine import calcfunction
from aiida import orm

@calcfunction
def add_hydrogens_to_structure(
    structure_data: orm.StructureData,
    potential_array: orm.ArrayData,
    do_supercell: orm.Bool, 
    equiv_peak_threshold: orm.Float,
    num_H: orm.Int
    ) -> dict:
    """Add hydrogen atoms to a structure based on its calculated potential."""

    def find_peaks(potential, do_supercell, equiv_peak_threshold):
        if do_supercell:
                  # We search for peaks in a 3x3x3 supercell to find also peaks close to a cell edge
            supercell_size = (3,3,3)
            supercell_potential = np.tile(potential, supercell_size) 
            peak_locations = peak_local_max(supercell_potential, min_distance=3, exclude_border=False)

            supercell_mask = (np.floor_divide(peak_locations, potential.shape) == np.array([1,1,1])).all(axis=1)
            peak_locations = peak_locations[supercell_mask, :]
           
            peak_locations_orig = np.remainder(peak_locations, potential.shape)

        else:
            peak_locations_orig = peak_local_max(potential, min_distance=3, exclude_border=False)

        peak_values_orig = potential[tuple(peak_locations_orig.T)]

        #To check if the maxima are sorted:
        assert np.allclose(peak_values_orig, sorted(peak_values_orig)[::-1]), "Peak values are not sorted!"
        num_orig_peaks = len(peak_values_orig)

        # Filter down only to largest equivalent peaks, within a threshold 
        
        equiv_peak_mask = peak_values_orig > peak_values_orig[0] * equiv_peak_threshold   
        peak_locations = peak_locations_orig[equiv_peak_mask]
        peak_values = peak_values_orig[equiv_peak_mask]  #equiv_peak_mask FUNZIONA PROPRIO DA MASCHERA

        return peak_locations, peak_values, num_orig_peaks, peak_locations_orig, peak_values_orig

    new_structure = structure_data.get_pymatgen()
    potential = potential_array.get_array('data')

    #Now I want to apply the above function to the data = look for maxima
    peak_locations, _, _, peak_locations_orig, peak_values_orig = find_peaks(
        potential, do_supercell.value, equiv_peak_threshold.value
    ) 

    for scaled_pos in np.divide(peak_locations, potential.shape):

        if int(new_structure.composition['H']) == num_H.value:
            break

        new_structure.append('H', scaled_pos, validate_proximity=True)

    all_peaks = orm.ArrayData()
    all_peaks.set_array('peak_values', peak_values_orig)
    all_peaks.set_array('peak_positions', np.divide(peak_locations_orig, potential.shape))

    return {
        'new_structure': orm.StructureData(pymatgen=new_structure),
        'all_peaks': all_peaks
    }
