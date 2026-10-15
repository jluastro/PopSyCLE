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

def get_star_system_pos_mag(hdf5_file, filt='ubv_I', ext_law='Damineli16', recalc=True):
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
        
    ext_law : str
        extinction law for calculating apparent magnitudes

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

def count_stars_hdf5(hdf5_file, filt='ubv_I', mag_threshold=21, ext_law='Damineli16'):
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
        
    ext_law : str
        extinction law for calculating apparent magnitudes

    Returns
    -------
    stars_above_threshold : float
        Number of stars brighter than mag threshold (1e-6).
    """
    df_all_stars = get_star_system_pos_mag(hdf5_file, filt=filt, ext_law=ext_law)

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

def events_for_popclass(h5_file, max_stars_per_bin=3e3, add_dl=False):
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
    dls = []
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
                dls.append(dists[use_lens])
    print(f'Drew {len(np.concatenate(thetaEs))} events total')
    if not add_dl:
        return np.concatenate(rem_ids), np.concatenate(thetaEs), np.concatenate(piEs), np.concatenate(tEs), np.concatenate(weights)
    else:
        return np.concatenate(rem_ids), np.concatenate(thetaEs), np.concatenate(piEs), np.concatenate(tEs), np.concatenate(weights), np.concatenate(dls)

def statistical_eventrate(h5_file, solid_angle, filt='ubv_I', mag_threshold=21
                          max_stars_per_bin=1e4, tE_range=(0,np.inf),
                          ext_law='Damineli16'):
    """
    Statistically estimate event rates from a
    PopSyCLE singles-only catalog.
    
    Parameters
    ----------
    hdf5_file : str
        Filename of an hdf5 file.
        
    solid_angle : float
        area of catalog on-sky in square degrees
        
    filt :  str
        filter to use for source mag limit
        
    mag_threshold : float
        source mag limit

    max_stars_per_bin : str
        Maxinum number of stars per bin to use for the
        calculation. Prevents excessive memory/computation use.
        Default is 1e4.
        
    tE_range : (float,float)
        minimum and maximum tE in days for events to include
        in the rate estimation
        
    ext_law : str
        extinction law for calculating apparent magnitudes

    Returns
    -------
    tau : float
        microlensing optical depth
        
    gamma_area : float
        events per year square degree
        
    gamma_star : float
        events per year per source star

    avg_tE : float
        average event timescale
        
    n_lens : int
        total number of stars in the catalog
        
    n_source : int
        total number of catalog stars that meet the mag threshold
    """
    
    # Load up file and basic details
    hf = h5py.File(h5_file, 'r')
    stars_per_bin = count_stars_all_per_bin(h5_file)
    n_stars_all = np.sum(stars_per_bin)
    n_sources_all = count_stars_hdf5(h5_file, filt=filt,
            mag_threshold=mag_threshold, ext_law=ext_law)
    
    # Set up sums
    n_events = 0
    n_stars_used = 0
    n_sources_used = 0
    thetaE_murel_sum = 0.0
    thetaE2_sum = 0.0
    scale_stars_used = np.maximum(1,int(np.floor(np.max(stars_per_bin)/max_stars_per_bin)))
    
    # Iterate over bins
    for k in list(hf.keys()):
        if (k.startswith('l')) & ('_' not in k):
            print('running', k)
            dat = hf[k]

            if dat.shape[0] > 0:
                # Downsample and grab needed data
                patch = dat[::scale_stars_used]
                dists = patch['rad']
                mul = patch['mu_lcosb']
                mub = patch['mu_b']
                masses = patch['mass']
                mags = synthetic.calc_app_mag(dists, patch[filt],
                            patch['exbv'], filt_dict[filt][ext_law])
                rem_id_catalog = patch['rem_id']
                idx = np.arange(len(masses))
                idx_s = np.where(mags<mag_theshold)[0]
                del patch

                # Set up source and lens indices for valid pairs
                src_idxs, lens_idxs = np.meshgrid(idx_s, idx)
                src_idxs, lens_idxs = src_idxs.ravel(), lens_idxs.ravel()
                dist_comp = (dists[src_idxs] > dists[lens_idxs]) #source further than lens
                use_srcs = src_idxs[dist_comp]
                use_lens = lens_idxs[dist_comp]

                # Microlensing math
                c, G, mSun, pctom = 299792458, 6.6743e-11, 1.98840987e+30, 3.08567758e+16
                theta_e = np.sqrt(4*G*mSun*masses[use_lens]*pi_rel/(1000*pctom*c**2)) #rad
                mu_rel = np.sqrt((mul[use_lens]-mul[use_srcs])**2 + (mub[use_lens]-mub[use_srcs])**2) * (1./1000.) * (1./365.25) * (1.0/60.0/60.0) * np.pi/180.0 # rad/day
                
                # Accounting
                n_events += len(use_lens)
                n_stars_used += len(idx)
                n_sources_used += len(idx_s)
                thetaE_murel_sum += np.sum(theta_e*murel)
                thetaE2_sum += np.sum(theta_e**2)
                
    print(f'Drew {n_events} events total')
    area_eff_lens = solid_angle * n_stars_used/n_stars_all
    area_eff_source = solid_angle * n_sources_used/n_sources_all
    # Optical depth
    tau = np.pi*thetaE2_sum / (n_sources_used * area_eff_lens)
    # Event rates
    gamma_area = 2*thetaE_murel_sum / area_eff_lens / (area_eff_source /
                    (np.pi/180)**2) * 365.25
    gamma_star = 2*thetaE_murel_sum / area_eff_lens / n_sources_used * 365.25
    # Average timescale
    avg_tE = thetaE2_sum / thetaE_murel_sum
    
    return tau, gamma_area, gamma_star, avg_tE, n_stars_all, n_sources_all


def _calc_blend_clusters(points, radius, star_idxs):
    """
    Helper function for calc_blends to run the KDTree clustering
    """
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

def _calc_blends_bin(star_dat, blend_rad, filters,
                primary_filter=None, ext_law='Damineli16'):
    """
    Helper function for calc_blends to process 1 bin
    """
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
    all_pts = np.transpose([delta_l_cosb, delta_b])
    star_idxs = star_dat['obj_id'].to_numpy()

    # Set up arrays with necessary data
    mags = star_dat[filters].to_numpy()
    mags_prim = star_dat[primary_filter].to_numpy()
    star_dat.loc[:,'plx'] = 1/star_dat['rad']
    astrom_vals = star_dat[['glon','glat','mu_lcosb',
                           'mu_b','plx','exbv']].to_numpy()

    # Run the cluster calculation on the sorted/cleaned points
    clusters0, prim_idxs = _calc_blend_clusters(all_pts, blend_rad, star_idxs)
    # Turn this into a masked array for easy mag + param sums
    clusters_arr = np.array(list(zip_longest(*clusters0, fillvalue=-1))).T
    clusters = np.ma.masked_values(clusters_arr, -1)

    # Compute the magnitude sums & weighted astrometry parameters
    mags = np.append(mags, [np.repeat(np.inf, len(filters))], axis=0)
    bmags = -2.5*np.log10(np.sum(10**(-0.4*mags[clusters,:]),axis=1))
    out = {}
    for i,f in enumerate(filters):
        out[f] = bmags[:,i]
    mags_prim = np.append(mags_prim, [np.inf])
    astrom_vals = np.append(astrom_vals, [np.repeat(0, 6)], axis=0)
    fluxes = 10**(-0.4*mags_prim[clusters])
    bvals = (np.sum(fluxes.T*astrom_vals[clusters,:].T, axis=1)
                / np.sum(fluxes,axis=1)).T
    for i, col in enumerate(['glon','glat','mu_lcosb',
                             'mu_b','plx', 'exbv']):
        out[col] = bvals[:,i]
    out['obj_id'] = prim_idxs
    return out

def calc_blends(hdf5_file, blend_rad, filters,
                primary_filter=None, ext_law='Damineli16',
                recalc=False, combine_bins=False):
    """
    Artificially blend the catalog based on a given radius &
    save the blend catalog to a file. Outputs are: summed magnitudes
    in selected filters, primary filter flux-weighted positions,
    proper motions, parallaxes, and reddening, and the obj_id in
    the original table of the brightest star in the blend.

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
    ext_law : str = 'Damineli16'
        extinction law
    recalc : boolean = False
        if True, recalculate and replace existing file
    combine_bins : boolean = False
        if True, save a single dataset instead of binned data
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
            blend_bin = _calc_blends_bin(pd.DataFrame(hf[dset_name][:]),
                blend_rad, filters,
                primary_filter=primary_filter, ext_law=ext_law)
            compound_dtype = synthetic._generate_compound_dtype(blend_bin)
            save_data = np.empty(len(blend_bin['plx']), dtype=compound_dtype)
            for colname in blend_bin:
                save_data[colname] = blend_bin[colname]
            if combine_bins:
                if 'data' not in outfile:
                    dataset = outfile.create_dataset('data', shape=(0,),
                                                chunks=(1e4,),
                                                maxshape=(None,),
                                                dtype=compound_dtype)
                else:
                    dataset = outfile['data']
                old_size = dataset.shape[0]
            else:
                if dset_name not in outfile:
                    dataset = outfile.create_dataset(dset_name, shape=(0,),
                                                chunks=(1e4,),
                                                maxshape=(None,),
                                                dtype=compound_dtype)
                else:
                    dataset = outfile[dset_name]
                old_size = 0
            new_size = old_size + len(blend_bin['plx'])
            dataset.resize((new_size, ))
            dataset[old_size:new_size] = save_data
    hf.close()
    outfile.close()
