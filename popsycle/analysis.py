import pandas as pd
import numpy as np
import h5py
import pylab as plt

from popsycle import utils
from popsycle import synthetic
import time
import os

def get_star_system_pos_mag(hdf5_file, filt='ubv_I', recalc=True):
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
    outfile = hdf5_file.replace('.h5', '_stars_posmag.pkl')

    if os.path.exists(outfile) and recalc == False:
        start_time = time.time()
        df_final = pd.read_pickle(outfile)
        stop_time = time.time()
        print(f'get_star_system_pos_mag: Run time = {stop_time - start_time} sec with recalc=False')
        return df_final

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

    list_of_df = []

    # Loop through the fields and aggregate the stars.
    print('Countint stars in patches: ', end="")
    for k in field_keys:
        print(f'{k}, ', end="")
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

            # Save to list of all data frames (to be concatenatted later)
            list_of_df.append(patch_df)

    df_final = pd.concat(list_of_df)
    del list_of_df

    stop_time = time.time()
    print()
    print(f'get_star_system_pos_mag: Run time = {stop_time - start_time} sec')

    df_final.to_pickle(outfile)

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

def count_stars_all(h5_file):
    """
    Finds the number of stars or systems in the field.

    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.

    Returns
    -------
    n_stars : int
        Number of stars.
    """
    hf = h5py.File(h5_file, 'r')
    n_stars = 0
    for k in list(hf.keys()):
        if '_' not in k:
            n_stars += hf[k].shape[0]
    return n_stars

def count_stars_all_per_bin(h5_file):
    """
    Finds the number of stars or systems per bin.

    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.

    Returns
    -------
    n_stars : list
        Number of stars per bin.
    """
    hf = h5py.File(h5_file, 'r')
    n_stars = []
    for k in list(hf.keys()):
        if '_' not in k:
            n_stars.append(hf[k].shape[0])
    return n_stars

def events_for_popclass(h5_file, max_stars_per_bin=3e3):
    """
    Draw microlensing events for random lens, source pairs from a 
    PopSyCLE singles-only catalog. 
    
    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.

    use_stars_per_bin : str
        Order of magnitude number of stars per bin to use for the 
        calculation. Prevents excessive memory/computation use.
        Default is 1e3.

    Returns
    -------
    thetaEs : np.array
        array of Einstein ring radii for events (mas)
    piEs : np.array
        array of microlensing parallaxes for events (unitless)
    tEs : np.array
        array of timescales for events (days)
    weights : np.array
        array of relative weights for events (mu_rel * thetaE)
    """
    hf = h5py.File(h5_file, 'r')
    thetaEs = []
    piEs = []
    tEs = []
    weights = []
    stars_per_bin = count_stars_all_per_bin(h5_file)
    scale_stars_used = np.maximum(1,int(np.floor(np.max(stars_per_bin)/max_stars_per_bin)))
    for k in list(hf.keys()):
        if '_' not in k:
            print('running', k)
            dat = hf[k]

            if dat.shape[0] > 0:
                patch = dat[::scale_stars_used]
                dists = patch['rad']
                mul = patch['mu_lcosb']
                mub = patch['mu_b']
                masses = patch['mass']
                idx = np.arange(len(masses))
                del patch

                src_idxs, lens_idxs = np.meshgrid(idx, idx)
                src_idxs, lens_idxs = src_idxs.ravel(), lens_idxs.ravel()

                dist_comp = (dists[src_idxs] > dists[lens_idxs]) #source further than lens
                use_srcs = src_idxs[dist_comp]
                use_lens = lens_idxs[dist_comp]

                # Microlensing math
                pi_rel = (1/dists[use_lens] - 1/dists[use_srcs])
                c, G, mSun, pctom = 299792458, 6.6743e-11, 1.98840987e+30, 3.08567758e+16
                theta_e = np.sqrt(4*G*mSun*masses[use_lens]*pi_rel/(1000*pctom*c**2)) * 180/np.pi * 60**2 * 1000
                pi_e = pi_rel / theta_e
                mu_rel = np.sqrt((mul[use_lens]-mul[use_srcs])**2 + (mub[use_lens]-mub[use_srcs])**2)
                t_e = theta_e/mu_rel * 365.25 # years -> days
                thetamu = theta_e*mu_rel
                thetaEs.append(theta_e)
                piEs.append(pi_e)
                tEs.append(t_e)
                weights.append(thetamu)
    print(f'Drew {len(np.concatenate(thetaEs))} events total')
    return np.concatenate(thetaEs), np.concatenate(piEs), np.concatenate(tEs), np.concatenate(weights)

