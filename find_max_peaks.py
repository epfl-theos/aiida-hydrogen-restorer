def find_max_peaks(filecube, do_supercell, equiv_peak_threshold, num_H):

    import time
    from pymatgen.io.common import VolumetricData
    from skimage.feature.peak import peak_local_max
    import numpy as np
    from ase.io import cube

    def find_peaks(potential, do_supercell, equiv_peak_threshold):
        import time
        from pymatgen.io.common import VolumetricData
        from skimage.feature.peak import peak_local_max
        import numpy as np
        from ase.io import cube
        
        if do_supercell:
               # We search for peaks in a 3x3x3 supercell to find also peaks close to a cell edge
            supercell_size = (3,3,3)
            supercell_potential = np.tile(potential, supercell_size) # to copy potential data inside the supercell

            peak_locations = peak_local_max(supercell_potential, min_distance=3, exclude_border=False)

            #Now I have the peaks in the 3x3x3 supercell but its necessary to rescale 
            # the data to have valid indexes for the original potential
            supercell_mask = (np.floor_divide(peak_locations, potential.shape) == np.array([1,1,1])).all(axis=1)
            peak_locations = peak_locations[supercell_mask, :]
           
            peak_locations_orig = np.remainder(peak_locations, potential.shape)

        else:
            peak_locations_orig = peak_local_max(potential, min_distance=3, exclude_border=False)

        #now we want the values of the peaks knowing their positions from peak_locations
        peak_values_orig = potential[tuple(peak_locations_orig.T)]

        #To check if the maxima are sorted:
        assert np.allclose(peak_values_orig, sorted(peak_values_orig)[::-1]), "Peak values are not sorted!"
        num_orig_peaks = len(peak_values_orig)

        # Filter down only to largest equivalent peaks, within a threshold 
        
        equiv_peak_mask = peak_values_orig > peak_values_orig[0] * equiv_peak_threshold   
        peak_locations = peak_locations_orig[equiv_peak_mask]
        peak_values = peak_values_orig[equiv_peak_mask]  #equiv_peak_mask FUNZIONA PROPRIO DA MASCHERA

        return peak_locations, peak_values, num_orig_peaks, peak_locations_orig, peak_values_orig
 
        

    t1 = time.time()

    # Read cube files (structure + grid data)

    cubefile = VolumetricData.from_cube('' + filecube)
    potential = cubefile.data['total']
    
    print(f"# Time to read cube file(1): {(time.time() - t1) * 1000:.2f} ms")
  
    cont = len(cubefile.structure.indices_from_symbol('H'))

    

    t1 = time.time()
    #Now I want to apply the above function to the data = look for maxima
    peak_locations, peak_values, num_orig_peaks, peak_locations_orig, peak_values_orig = find_peaks(potential, do_supercell, equiv_peak_threshold) 

    print(f"# Time to find peak locations: {(time.time() - t1) * 1000:.2f} ms")
    
    
    c = 0 #setting the counter of NEW H to zero

    #Converting to scaled coordinates adding them to the structure
         
    for peak_location, peak_value, scaled_pos in zip(
        peak_locations, peak_values, np.divide(peak_locations, potential.shape)):

        if cont + c < num_H:
            
            cubefile.structure.append('H', scaled_pos, validate_proximity=True)
            c = c + 1 

            print(f"Peak found at scaled coordinates {scaled_pos}")
            print(f"    -> Array position: {peak_location}")
            print(f"    -> peak value: {peak_value} ({(peak_value / peak_values[0]) * 100:.4f}% of max peak)")
         
        if cont + c == num_H:
            print('WARNING! Reached max number of Hydrogens.') #Maybe it's better to distinguish when the number of maxima are cut because they were more
                                                               # than max possible number of H AND when the number maxima selected by the threshold was the same as
                                                               # the number of remaining H to put in the cell.
            break
            
    print(f"Keeping only {c} peaks out of {num_orig_peaks}")

    new_H = cont + c  #number of H in the 'new' structure
 
    allpos_H = [peak_values_orig, np.divide(peak_locations_orig, potential.shape)]  #all maxima found, without considering the threshold


    return cubefile.structure, new_H, allpos_H






        
