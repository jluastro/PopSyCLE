import pandas as pd
import numpy as np
import h5py
import pylab as plt
from scipy.spatial import KDTree
from popsycle import utils, synthetic, phot_utils
import time
import os
import pdb
import warnings
from itertools import zip_longest

filt_dict = phot_utils.make_filt_dict()

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
    print('Counting stars in patches: ', end="")
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
        if (k.startswith('l')) & ('_' not in k):
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
        if (k.startswith('l')) & ('_' not in k):
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

    max_stars_per_bin : str
        Maxinum number of stars per bin to use for the 
        calculation. Prevents excessive memory/computation use.
        Default is 3e3.

    Returns
    -------
    rem_ids : np.array
        array of lens types
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
    rem_ids = []
    thetaEs = []
    piEs = []
    tEs = []
    weights = []
    stars_per_bin = count_stars_all_per_bin(h5_file)
    scale_stars_used = np.maximum(1,int(np.floor(np.max(stars_per_bin)/max_stars_per_bin)))
    for k in list(hf.keys()):
        if (k.startswith('l')) & ('_' not in k):
            print('running', k)
            dat = hf[k]

            if dat.shape[0] > 0:
                patch = dat[::scale_stars_used]
                dists = patch['rad']
                mul = patch['mu_lcosb']
                mub = patch['mu_b']
                masses = patch['mass']
                rem_id_catalog = patch['rem_id']
                idx = np.arange(len(masses))
                del patch

                src_idxs, lens_idxs = np.meshgrid(idx, idx)
                src_idxs, lens_idxs = src_idxs.ravel(), lens_idxs.ravel()

                dist_comp = (dists[src_idxs] > dists[lens_idxs]) #source further than lens
                use_srcs = src_idxs[dist_comp]
                use_lens = lens_idxs[dist_comp]
                rem_id = rem_id_catalog[use_lens]

                # Microlensing math
                pi_rel = (1/dists[use_lens] - 1/dists[use_srcs])
                c, G, mSun, pctom = 299792458, 6.6743e-11, 1.98840987e+30, 3.08567758e+16
                theta_e = np.sqrt(4*G*mSun*masses[use_lens]*pi_rel/(1000*pctom*c**2)) * 180/np.pi * 60**2 * 1000
                pi_e = pi_rel / theta_e
                mu_rel = np.sqrt((mul[use_lens]-mul[use_srcs])**2 + (mub[use_lens]-mub[use_srcs])**2)
                t_e = theta_e/mu_rel * 365.25 # years -> days
                thetamu = theta_e*mu_rel
                rem_ids.append(rem_id)
                thetaEs.append(theta_e)
                piEs.append(pi_e)
                tEs.append(t_e)
                weights.append(thetamu)
    print(f'Drew {len(np.concatenate(thetaEs))} events total')
    return np.concatenate(rem_ids), np.concatenate(thetaEs), np.concatenate(piEs), np.concatenate(tEs), np.concatenate(weights)


# For the blend calculation, sum the mags
def blend_magsum(magss):
    return -2.5*np.log10(np.sum(10**(-0.4*magss),axis=1))

# For the blend calculation, flux-weight average the astrometric values
def blend_mag_weighted_sum(mags, vals):
    fluxes = 10**(-0.4*mags)
    sum_vals = (np.sum(fluxes.T*vals.T, axis=1)/np.sum(fluxes,axis=1)).T
    return sum_vals

# For the blend calculation, use a KDTree to assign stars to blend "clusters"
def calc_blend_clusters(points, radius, star_idxs):
    # Set up KDTree and array to note whether each point has been assigned a cluster
    tree = KDTree(points)
    n_points = len(points)
    assigned = np.zeros(n_points, dtype=bool)
    clusters = []
    prim_idxs = []
    for i in range(n_points):
        # Find unassigned neighbors and save cluster
        if not assigned[i]:
            neighbor_indices = tree.query_ball_point(points[i], r=radius)
            current_cluster = [idx for idx in neighbor_indices if not assigned[idx]]
            assigned[current_cluster] = True
            clusters.append(current_cluster)
            prim_idxs.append(star_idxs[i])
    return clusters, np.array(prim_idxs)

"""
Take a resolved star catalog and blend it, summing magnitudes
and calculating flux-weighted positions, proper
motions, and parallaxes.
"""
def calc_blends_bin(star_dat, blend_rad, filters,
                primary_filter=None, ext_law='Damineli16'):
    # Clean up table and sort by magnitude
    if primary_filter is None:
        warnings.warn(f"No primary_filter provided. Using {filters[0]}"+
            " to sort stars and flux-weight additional parameters.")
        primary_filter = filters[0]
    for filt in filters:
        star_dat.loc[:,filt] = synthetic.calc_app_mag(star_dat['rad'].to_numpy(),
                star_dat[filt].to_numpy(), star_dat['exbv'].to_numpy(),
                filt_dict[filt][ext_law])
    star_dat = star_dat[~np.isnan(star_dat[primary_filter])].copy()
    star_dat.sort_values(by=primary_filter, inplace=True)
    star_dat.reset_index(drop=True,inplace=True)

    # Get coordinates into simplified delta_l_cosb and delta_b frame
    l_mean = np.mean(star_dat['glon']); b_mean = np.mean(star_dat['glat'])
    all_l, all_b = star_dat['glon'].to_numpy(), star_dat['glat'].to_numpy()
    delta_l_cosb = ((all_l-360*(all_l>180))-l_mean)*np.cos(all_b*np.pi/180)
    delta_b = all_b-b_mean
    all_pts = np.transpose([all_l,all_b])
    star_idxs = star_dat['obj_id'].to_numpy()

    # Set up arrays with necessary data
    mags = star_dat[filters].to_numpy()
    mags_prim = star_dat[primary_filter].to_numpy()
    star_dat.loc[:,'plx'] = 1/star_dat['rad']
    astrom_vals = star_dat[['glon','glat','mu_lcosb',
                           'mu_b','plx','exbv']].to_numpy()
    blend_vals = []

    # Run the cluster calculation on the sorted/cleaned points
    clusters0, prim_idxs = calc_blend_clusters(all_pts, blend_rad, star_idxs)
    # Turn this into a masked array for easy mag + param sums
    clusters_arr = np.array(list(zip_longest(*clusters0, fillvalue=-1))).T
    clusters = np.ma.masked_values(clusters_arr, -1)

    # Compute the magnitude sums & other values if desired
    mags = np.append(mags, [np.repeat(np.inf, len(filters))], axis=0)
    bmags = blend_magsum(mags[clusters,:])
    out = {}
    for i,f in enumerate(filters):
        out[f] = bmags[:,i]
    mags_prim = np.append(mags_prim, [np.inf])
    astrom_vals = np.append(astrom_vals, [np.repeat(0, 6)], axis=0)
    bvals = blend_mag_weighted_sum(mags_prim[clusters], astrom_vals[clusters,:])
    for i, col in enumerate(['glon','glat','mu_lcosb',
                             'mu_b','plx', 'exbv']):
        out[col] = bvals[:,i]
    out['obj_id'] = prim_idxs
    print(out.keys())
    return out

def calc_blends(hdf5_file, blend_rad, filters,
                primary_filter=None, recalc=False):
    """
    Artificially blend the catalog based on a given radius &
    save the blend catalog to a file.

    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.
    blend_radius : float
        blend radius in degrees
    filters : list of str
        column headers for photometry to blend
    primary_filter : str = None
        filter to use for astrometric parameter flux-weighting.
        if None, use the first filter in the list.
    recalc : boolean = False
        if True, recalculate and replace existing file
    """
    outfile_name = hdf5_file.replace('.h5', '_blended_catalog.h5')

    if os.path.exists(outfile_name) and recalc == False:
        print(f'{outfile_name} exists. Set recalc=True to recalculate.')
        return
    
    hf = h5py.File(hdf5_file, 'r')
    outfile = h5py.File(outfile_name, 'a')

    for dset_name in list(hf.keys()):
        if (dset_name.startswith('l')) & ('_' not in dset_name):
            print(dset_name)
            blend_bin = calc_blends_bin(pd.DataFrame(hf[dset_name][:]),
                blend_rad, filters,
                primary_filter=primary_filter)
            compound_dtype = synthetic._generate_compound_dtype(blend_bin)
            save_data = np.empty(len(blend_bin['plx']), dtype=compound_dtype)
            for colname in blend_bin:
                save_data[colname] = blend_bin[colname]
            if dset_name not in outfile:
                dataset = outfile.create_dataset(dset_name, shape=(0,),
                                            chunks=(1e4,),
                                            maxshape=(None,),
                                            dtype=compound_dtype)
            else:
                dataset = outfile[dset_name]
            new_size = len(blend_bin['plx'])
            dataset.resize((new_size, ))
            dataset[:] = save_data
    hf.close()
    outfile.close()
