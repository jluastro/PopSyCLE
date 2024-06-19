import pandas as pd
import numpy as np
import h5py
import pylab as plt

import utils
from popsycle import synthetic
import time

def get_star_system_pos_mag(hdf5_file, filt='ubv_I'):
    """
    Return a table with lists of star systems and their RA, Dec, z position,
    and system apparent magnitude. This is useful for making stellar density maps,
    computing microlens event occurrence rates, etc. Columns will also be
    returned containing the

    Parameters
    ----------
    hdf5_file : str
        Name of the H5 file output from perform_pop_syn.
    filt : str
        filter to use to calculate the apparent magnitude

    Returns
    -------
    df_final : pandas dataframe
        Columns include:
            ['obj_id', 'exbv', 'glat', 'glon', 'rad', 'isMultiple',
            'N_companions', 'rem_id', 'm_ubv_I_app']

        where the magnitude column is the apparent system magnitude
        in the designated filter.
    """
    # Load up the H5 file with star systems.
    hf = h5py.File(hdf5_file, 'r')

    # Trim keys down to valid patch keys (e.g. "l0b0").
    orig_keys = hf.keys()   # Patch keys.
    field_keys = []
    for key in list(orig_keys)[:-2]:
        if key.startswith('l'):
            field_keys.append(key)

    ext_law = 'Damineli16'
    start_time = time.time()
    df_final = None

    # Loop through the fields and aggregate the stars.
    for k in field_keys:
        patch = np.array(hf[k])

        if len(patch) > 0:
            # Make a Pandas data frame. Faster to work with and index against with companions.
            patch_df = pd.DataFrame(data=patch, columns=np.dtype(patch[0]).names)
            patch_df.set_index(['obj_id'])

            # Memory management
            del patch
            patch_df.drop(columns=['px', 'py', 'pz', 'vx', 'vy', 'vz',
                                   'vr', 'mu_b', 'mu_lcosb',
                                   'zams_mass', 'mass', 'systemMass',
                                   'age', 'popid', 'mbol', 'grav', 'teff', 'feh', 'mbol'],
                          inplace=True)

            # Make flux column.
            patch_df['m_' + filt + '_app'] = synthetic.calc_app_mag(patch_df['rad'],
                                                                    patch_df[filt],
                                                                    patch_df['exbv'],
                                                                    synthetic.filt_dict[filt][ext_law])

            # More memory management. Drop absolute mag columns
            patch_df.drop(patch_df.filter(regex='^ubv').columns, axis=1, inplace=True)

            # Save to final stacked dataframe (all stars)
            if df_final is None:
                df_final = patch_df
            else:
                df_final = pd.concat([df_final, patch_df])

    del hf

    stop_time = time.time()
    print(f'Run time = {stop_time - start_time} sec')

    return df_final

def count_stars_hdf5(hdf5_file, filt='ubv_I', mag_threshold=21):
    """
    Finds the number of stars in the field brighter than a certain mag.
    Assumes binary/multiple stars are blended.

    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.

    filt : str
        Ubv filter in PopSyCLE ubv_(U, B, V, I, R, J, H, K).
        Default is ubv_I.

    mag_threshold: float
        Magnitude below which we count the number of stars.

    Returns
    -------
    stars_above_threshold : float
        Number of stars brighter than mag threshold (1e-6).
    """
    df_all_stars = get_star_system_pos_mag(hdf5_file, filt=filt)

    # Get good stars above our magnitude threshold.
    gdx = np.where(df_all_stars['m_' + filt + '_app'] < mag_threshold)[0]
    N_stars = len(gdx)

    return N_stars
