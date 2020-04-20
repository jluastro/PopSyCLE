#! /usr/bin/env python
"""
converter.py
Functions for converting PopSyCLE data formats.
"""

import h5py
import numpy as np
from popsycle.synthetic import _generate_comp_dtype

def convert_h5_array_dtype_to_compound_dtype(hdf5_file):
    """
    Converts h5 files with the array datatype into the compound datatype

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file with the array datatype.

    Output
    ------
    <hdf5_file>.compound_dtype.fits
        Name of the resulting HDF5 file with the compound datatype.

    """
    # Load in the original col_idx
    col_idx = {'zams_mass': 0, 'rem_id': 1, 'mass': 2,
               'px': 3, 'py': 4, 'pz': 5,
               'vx': 6, 'vy': 7, 'vz': 8,
               'rad': 9, 'glat': 10, 'glon': 11,
               'vr': 12, 'mu_b': 13, 'mu_lcosb': 14,
               'age': 15, 'popid': 16, 'ubv_K': 17,
               'ubv_I': 18, 'exbv': 19, 'obj_id': 20,
               'ubv_J': 21, 'ubv_U': 22, 'ubv_R': 23,
               'ubv_B': 24, 'ubv_H': 25, 'ubv_V': 26}

    # Create a reverse col_idx dictionary where the column number
    # points to the column name
    reverse_col_idx = {}
    for key, val in col_idx.items():
        reverse_col_idx[val] = key

    # Build a compound datatype that assigns a variable type for each key
    comp_dtype = _generate_comp_dtype(col_idx)

    # Load in the old file and prepare to write to the new file
    hdf5_file_new = hdf5_file.replace('.h5', '.compound_dtype.h5')
    f_in = h5py.File(hdf5_file, 'r')
    f_out = h5py.File(hdf5_file_new, 'w')

    # Looping over all of the datasets
    for key in f_in:
        # Simply copy over the 'lat_bin_edges' and 'long_bin_edges' datasets
        if 'bin_edges' in key:
            f_out[key] = f_in[key][:]
            continue

        # Load in the data
        data_in = f_in[key][:]
        # Prepare an empty data_out array with the compound datatype
        data_out = np.empty(data_in.shape[1], dtype=comp_dtype)
        # Load up the data_out array with the columns
        # given the names in the reverse_col_idx
        for i in range(data_in.shape[0]):
            data_out[reverse_col_idx[i]] = data_in[i, :]

        # Write that data out to the new file and free up the data_in
        f_out[key] = data_out
        del data_in

    f_in.close()
    f_out.close()

    print('%s generated' % hdf5_file_new)
