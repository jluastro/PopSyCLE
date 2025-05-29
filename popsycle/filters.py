#! /usr/bin/env python
"""
filters.py
Functions (and their associated functions) for additional photometric systems.
"""
import inspect
import os
import numpy as np
from spisea import evolution, atmospheres, reddening, ifmr, synthetic
from spisea.imf import multiplicity, imf
from scipy.interpolate import griddata
from scipy.spatial.ckdtree import cKDTree


def generate_ubv_to_ztf_grid(iso_dir, filter_name):
    """
    Creates the 2D transformational matrix `ubv_to_ztf-r_grid.npz` and
    `ubv_to_ztf-g_grid.npz` necessary for generating ztf-g and ztf-r
    magnitudes from the UBV filters.

    2D transformational matrix is valid for values of:
        0 < ubv_V - ubv_R < 6
        0 < ubv_B - ubv_V < 6

    ubv-to-ztf-g
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_g - ubv_V

    ubv-to-ztf-r
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_r - ubv_R

    ubv-to-ztf-i
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_i - ubv_I

    Parameters
    ----------
    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. Must be either 'g' or 'r'.

    """

    # Define isochrone parameters for calculating absolute magnitudes
    logAge = np.log10(8 * 10 ** 9)  # Age in log(years)
    dist = 10  # distance in parsec
    metallicity = 0  # Metallicity in [M/H]

    # Define evolution/atmosphere models and extinction law
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawDamineli16()

    # Also specify filters for synthetic photometry
    filt_list = ['ztf,r', 'ztf,g', 'ztf,i',
                 'ubv,B', 'ubv,V', 'ubv,R', 'ubv,I']
    filt_list_reformat = ['m_%s' % f.replace(',', '_') for f in filt_list]

    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Kroupa+01
    # Define boundaries of each mass segement
    massLimits = np.array([0.08, 0.5, 1, 120])
    # Power law slope associated with each mass segment
    powers = np.array([-1.3, -2.3, -2.3])
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)

    # Define total cluster mass
    mass = 10 ** 5.

    ubv_b = np.array([])
    ubv_v = np.array([])
    ubv_r = np.array([])
    ubv_i = np.array([])
    ztf_g = np.array([])
    ztf_r = np.array([])
    ztf_i = np.array([])

    # Create photometry for a range of extinctions
    for AKs in np.arange(0, 1.1, .1):
        my_iso = synthetic.IsochronePhot(logAge, AKs, dist,
                                         metallicity=metallicity,
                                         evo_model=evo_model,
                                         atm_func=atm_func,
                                         red_law=red_law,
                                         filters=filt_list,
                                         iso_dir=iso_dir)
        # Check that the isochrone has all of the filters in filt_list
        # If not, force recreating the isochrone with recomp=True
        iso_filters = [f for f in my_iso.points.colnames if 'm_' in f]
        if len(set(filt_list_reformat) - set(iso_filters)) > 0:
            my_iso = synthetic.IsochronePhot(logAge, AKs, dist,
                                             metallicity=metallicity,
                                             evo_model=evo_model,
                                             atm_func=atm_func,
                                             red_law=red_law,
                                             filters=filt_list,
                                             iso_dir=iso_dir,
                                             recomp=True)
        # Make cluster object
        cluster = synthetic.ResolvedCluster(my_iso, my_imf, mass,
                                            ifmr=None)
        clust = cluster.star_systems
        cond = ~np.isnan(clust['m_ubv_V'])
        clust = clust[cond]
        clust_cond = np.random.choice(np.arange(len(clust)),
                                      size=10000, replace=False)

        ubv_b = np.append(ubv_b, clust['m_ubv_B'][clust_cond])
        ubv_v = np.append(ubv_v, clust['m_ubv_V'][clust_cond])
        ubv_r = np.append(ubv_r, clust['m_ubv_R'][clust_cond])
        ubv_i = np.append(ubv_i, clust['m_ubv_I'][clust_cond])
        ztf_g = np.append(ztf_g, clust['m_ztf_g'][clust_cond])
        ztf_r = np.append(ztf_r, clust['m_ztf_r'][clust_cond])
        ztf_i = np.append(ztf_i, clust['m_ztf_i'][clust_cond])

    # Given the filter name, define a difference in magnitude to be fit for
    if filter_name == 'g':
        delta_m = ztf_g - ubv_v
    elif filter_name == 'r':
        delta_m = ztf_r - ubv_r
    elif filter_name == 'i':
        delta_m = ztf_i - ubv_i


    # Colors in both x and y direction go from 0 to 6 magnitudes
    # x_grid_arr: ubv_v - ubv_r
    # y_grid_arr: ubv_b - ubv_v
    x_grid_arr = np.linspace(0, 6, 1000)
    y_grid_arr = np.linspace(0, 6, 1000)

    # Create a grid of values on x_grid_arr and y_grid_arr
    # with linear algorithm
    ubv_to_ztf_grid = griddata((ubv_v - ubv_r, ubv_b - ubv_v),
                               delta_m,
                               (x_grid_arr[None, :], y_grid_arr[:, None]),
                               method='linear')

    # Resample this grid with both the liner and nearest algorithms onto a
    # finer grid. This allows for the 'nearest' method to
    # create fewer artifacts
    xx, yy = np.meshgrid(x_grid_arr, y_grid_arr)
    xx, yy = xx.flatten(), yy.flatten()

    cond = ~np.isnan(ubv_to_ztf_grid.flatten())
    ubv_to_ztf_grid_filled = griddata((xx[cond], yy[cond]),
                                      ubv_to_ztf_grid.flatten()[cond],
                                      (x_grid_arr[None, :],
                                       y_grid_arr[:, None]),
                                      method='linear')

    ubv_to_ztf_grid_nearest = griddata((xx[cond], yy[cond]),
                                       ubv_to_ztf_grid.flatten()[cond],
                                       (x_grid_arr[None, :],
                                        y_grid_arr[:, None]),
                                       method='nearest')

    # Place values into final grid from linear algorthm, and from the
    # nearest algorithm where the linear algorithm could not find a solution
    cond = np.isnan(ubv_to_ztf_grid_filled)
    ubv_to_ztf_grid_final = np.zeros_like(ubv_to_ztf_grid_filled)
    ubv_to_ztf_grid_final[cond] = ubv_to_ztf_grid_nearest[cond]
    ubv_to_ztf_grid_final[~cond] = ubv_to_ztf_grid_filled[~cond]

    # Save the data
    grid_arr = np.squeeze(np.dstack([xx, yy]), axis=0)
    data_dir = '%s/data' % os.path.dirname(inspect.getfile(generate_ubv_to_ztf_grid))
    ubv_to_ztf_filename = '%s/ubv_to_ztf-%s_grid.npz' % (data_dir, filter_name)
    np.savez(ubv_to_ztf_filename,
             ubv_to_ztf_grid=ubv_to_ztf_grid_final.astype(np.float32),
             kdtree_grid=grid_arr.astype(np.float32))


def load_ubv_to_ztf_grid(filter_name):
    """
    Loads the 2D transformational matrix `ubv_to_ztf-r_grid.npz` and
    `ubv_to_ztf-g_grid.npz` necessary for generating ztf-g and ztf-r
    magnitudes from the UBV filters, as well as the kdtree of those values.

    2D transformational matrix is valid for values of:
        0 < ubv_V - ubv_R < 6
        0 < ubv_B - ubv_V < 6

    ubv-to-ztf-g
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_g - ubv_V

    ubv-to-ztf-r
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_r - ubv_R

    ubv-to-ztf-i
        x-axis : ubv_V - ubv_R
        y-axis : ubv_B - ubv_V
        z-axis : ztf_i - ubv_I

    Parameters
    ----------
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. Must be either 'g' or 'r' or 'i'.

    Returns
    -------
    ubv_to_ztf_grid : 2D numpy array
        2D grid array of UBV colors with each cell containing the difference
        between a ztf filter and a ubv filter

    kdtree : cKDTree
        kdtree containing the grid of colors on the x-axis and y-axis

    """
    # Check for correct filter
    if filter_name not in ['g', 'r', 'i']:
        raise Exception("filter_name must be in: ['g', 'r', 'i']")

    # Load the ubv_to_ztf_grid from the file
    data_dir = '%s/data' % os.path.dirname(inspect.getfile(load_ubv_to_ztf_grid))
    ubv_to_ztf_filename = '%s/ubv_to_ztf-%s_grid.npz' % (data_dir, filter_name)
    ubv_to_ztf_grid_file = np.load(ubv_to_ztf_filename)

    # Generate a kdtree at the locations of all of the grid points
    ubv_to_ztf_grid = ubv_to_ztf_grid_file['ubv_to_ztf_grid']
    kdtree = cKDTree(ubv_to_ztf_grid_file['kdtree_grid'])

    return ubv_to_ztf_grid, kdtree


def transform_ubv_to_ztf(filter_name, ubv_B, ubv_V, ubv_R, ubv_I=None):
    """
    Converts ubv filters into ztf filters.

    Function is valid for values of:
        0 < ubv_V - ubv_R < 6
        0 < ubv_B - ubv_V < 6

    Parameters
    ----------
    filter_name : str
        ZTF filter name of converted photometry
        Can be 'g', 'r', or 'i'
        If converting to ztf_i, then ubv_I must be provided

    ubv_B : array of floats
        ubv_B photometry of galaxia / SPISEA sources

    ubv_V : array of floats
        ubv_V photometry of galaxia / SPISEA sources

    ubv_R : array of floats
        ubv_R photometry of galaxia / SPISEA sources

    ubv_I : array of floats, optional
        ubv_I photometry of galaxia / SPISEA sources,
        required for converting to ztf_i

    Returns
    -------
    ztf_mag : array of floats
        ztf photometry of galaxia / SPISEA sources in `filter_name`

    """
    # Check for correct filter
    if filter_name not in ['g', 'r', 'i']:
        raise Exception("filter_name must be in: ['g', 'r', 'i']")

    if filter_name == 'i' and ubv_I is None:
        raise Exception('ubv_I must be provided to convert to ztf_i')

    # Convert the ubv photometry into the right format
    x_data = ubv_V - ubv_R
    y_data = ubv_B - ubv_V
    data = np.squeeze(np.dstack([x_data, y_data]), axis=0)

    # Only query on data that is luminous
    cond_lum = ~np.isnan(data).any(axis=1)

    # Start with an empty array of nans
    ztf_diff = np.full(len(ubv_B), np.nan)

    # Find locations on the grid where x_data and y_data are located.
    # Put those values into the ztf_diff arrays
    ubv_to_ztf_grid, kdtree = load_ubv_to_ztf_grid(filter_name)
    _, indexes = kdtree.query(data[cond_lum])
    ztf_diff[cond_lum] = ubv_to_ztf_grid.flatten()[indexes]

    # Convert to ztf_g, ztf_r or ztf_i
    if filter_name == 'g':
        ztf_mag = ubv_V + ztf_diff
    elif filter_name == 'r':
        ztf_mag = ubv_R + ztf_diff
    elif filter_name == 'i':
        ztf_mag = ubv_I + ztf_diff

    return ztf_mag


def ztf_mag_vega_to_AB(ztf_mag_vega, filter_name):
    """
    Converts vega magnitudes into AB magnitudes for ztf filters.
    Extrapolated from http://astroweb.case.edu/ssm/ASTR620/alternateabsmag.html
    using the effective wavelengths from
    http://svo2.cab.inta-csic.es/svo/theory/fps3/

    Parameters
    ----------
    ztf_mag_vega : float, array of floats
        ztf photometry of galaxia / SPISEA sources in vega system

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. Must be either 'g' or 'r' or 'i'.

    Returns
    -------
    ztf_mag_AB : float, array of floats
        ztf photometry of galaxia / SPISEA sources in AB system

    """
    if filter_name == 'g':
        ztf_mag_AB = ztf_mag_vega - 0.07
    elif filter_name == 'r':
        ztf_mag_AB = ztf_mag_vega + 0.19
    elif filter_name == 'i':
        ztf_mag_AB = ztf_mag_vega + 0.43
    else:
        raise Exception('filter_name must be either g or r or i')

    return ztf_mag_AB


def ztf_mag_AB_to_vega(ztf_mag_AB, filter_name):
    """
    Converts AB magnitudes into vega magnitudes for ztf filters.
    Extrapolated from http://astroweb.case.edu/ssm/ASTR620/alternateabsmag.html
    using the effective wavelengths from
    http://svo2.cab.inta-csic.es/svo/theory/fps3/

    Parameters
    ----------
    ztf_mag_AB : float, array of floats
        ztf photometry of galaxia / SPISEA sources in AB system

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. Must be either 'g' or 'r' or 'i'.

    Returns
    -------
    ztf_mag_vega : float, array of floats
        ztf photometry of galaxia / SPISEA sources in vega system

    """
    if filter_name == 'g':
        ztf_mag_vega = ztf_mag_AB + 0.07
    elif filter_name == 'r':
        ztf_mag_vega = ztf_mag_AB - 0.19
    elif filter_name == 'i':
        ztf_mag_vega = ztf_mag_AB - 0.43
    else:
        raise Exception('filter_name must be either g or r or i')

    return ztf_mag_vega

def generate_ubv_to_rubin_grid(iso_dir,filter_name):

    logAge = np.log10(8 * 10 ** 9)  # Age in log(years)
    dist = 10  # distance in parsec
    metallicity = 0  # Metallicity in [M/H]

    # Define evolution/atmosphere models and extinction law
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawDamineli16()

    # Also specify filters for synthetic photometry
    filt_list = ['rubin,u', 'rubin,g', 'rubin,r', 'rubin,i', 'rubin,z', 'rubin,y',
                 'ubv,U', 'ubv,B', 'ubv,V', 'ubv,R', 'ubv,I','2mass,J']
    filt_list_reformat = ['m_%s' % f.replace(',', '_') for f in filt_list]

    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Kroupa+01
    # Define boundaries of each mass segement
    massLimits = np.array([0.08, 0.5, 1, 120])
    # Power law slope associated with each mass segment
    powers = np.array([-1.3, -2.3, -2.3])
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)

    # Define total cluster mass
    mass = 10 ** 5.

    ubv_u = np.array([])
    ubv_b = np.array([])
    ubv_v = np.array([])
    ubv_r = np.array([])
    ubv_i = np.array([])
    twomass_j = np.array([])
    rubin_u = np.array([])
    rubin_g = np.array([])
    rubin_r = np.array([])
    rubin_i = np.array([])
    rubin_z = np.array([])
    rubin_y = np.array([])

    # Create photometry for a range of extinctions
    for AKs in np.arange(0, 1.1, .1):
        my_iso = synthetic.IsochronePhot(logAge, AKs, dist,
                                         metallicity=metallicity,
                                         evo_model=evo_model,
                                         atm_func=atm_func,
                                         red_law=red_law,
                                         filters=filt_list,
                                         iso_dir=iso_dir)
        # Check that the isochrone has all of the filters in filt_list
        # If not, force recreating the isochrone with recomp=True
        iso_filters = [f for f in my_iso.points.colnames if 'm_' in f]
        if len(set(filt_list_reformat) - set(iso_filters)) > 0:
            my_iso = synthetic.IsochronePhot(logAge, AKs, dist,
                                             metallicity=metallicity,
                                             evo_model=evo_model,
                                             atm_func=atm_func,
                                             red_law=red_law,
                                             filters=filt_list,
                                             iso_dir=iso_dir,
                                             recomp=True)
        # Make cluster object
        cluster = synthetic.ResolvedCluster(my_iso, my_imf, mass,
                                            ifmr=None)
        clust = cluster.star_systems
        cond = ~np.isnan(clust['m_ubv_V'])
        clust = clust[cond]
        clust_cond = np.random.choice(np.arange(len(clust)),
                                      size=10000, replace=False)

        ubv_u = np.append(ubv_u, clust['m_ubv_U'][clust_cond])
        ubv_b = np.append(ubv_b, clust['m_ubv_B'][clust_cond])
        ubv_v = np.append(ubv_v, clust['m_ubv_V'][clust_cond])
        ubv_r = np.append(ubv_r, clust['m_ubv_R'][clust_cond])
        ubv_i = np.append(ubv_i, clust['m_ubv_I'][clust_cond])
        twomass_j = np.append(twomass_j, clust['m_2mass_J'][clust_cond])
        rubin_u = np.append(rubin_u, clust['m_rubin_u'][clust_cond])
        rubin_g = np.append(rubin_g, clust['m_rubin_g'][clust_cond])
        rubin_r = np.append(rubin_r, clust['m_rubin_r'][clust_cond])
        rubin_i = np.append(rubin_i, clust['m_rubin_i'][clust_cond])
        rubin_z = np.append(rubin_z, clust['m_rubin_z'][clust_cond])
        rubin_y = np.append(rubin_y, clust['m_rubin_y'][clust_cond])

    # Given the filter name, define a difference in magnitude to be fit for
    if filter_name == 'u':
        delta_m = rubin_u - ubv_u
    elif filter_name == 'g':
        delta_m = rubin_g - ubv_v
    elif filter_name == 'r':
        delta_m = rubin_r - ubv_r
    elif filter_name == 'i':
        delta_m = rubin_i - ubv_i
    elif filter_name == 'z':
        delta_m = rubin_z - ubv_i
    elif filter_name == 'y':
        delta_m = rubin_y - twomass_j


    # Colors in both x and y direction go from 0 to 6 magnitudes
    # x_grid_arr: ubv_v - ubv_r
    # y_grid_arr: ubv_b - ubv_v
    x_grid_arr = np.linspace(0, 6, 1000)
    y_grid_arr = np.linspace(0, 6, 1000)

    # Create a grid of values on x_grid_arr and y_grid_arr
    # with linear algorithm
    if filter_name == 'u':
        x,y = (ubv_b - ubv_v, ubv_u - ubv_b)
    elif filter_name == 'y':
        x,y = (ubv_i - twomass_j, ubv_r - ubv_i)
    else:
        x,y = (ubv_v - ubv_r, ubv_b - ubv_v)
        
    ubv_to_rubin_grid = griddata((x,y), delta_m,
                               (x_grid_arr[None, :], y_grid_arr[:, None]),
                               method='linear')

    # Resample this grid with both the liner and nearest algorithms onto a
    # finer grid. This allows for the 'nearest' method to
    # create fewer artifacts
    xx, yy = np.meshgrid(x_grid_arr, y_grid_arr)
    xx, yy = xx.flatten(), yy.flatten()

    cond = ~np.isnan(ubv_to_rubin_grid.flatten())
    ubv_to_rubin_grid_filled = griddata((xx[cond], yy[cond]),
                                      ubv_to_rubin_grid.flatten()[cond],
                                      (x_grid_arr[None, :],
                                       y_grid_arr[:, None]),
                                      method='linear')

    ubv_to_rubin_grid_nearest = griddata((xx[cond], yy[cond]),
                                       ubv_to_rubin_grid.flatten()[cond],
                                       (x_grid_arr[None, :],
                                        y_grid_arr[:, None]),
                                       method='nearest')

    # Place values into final grid from linear algorthm, and from the
    # nearest algorithm where the linear algorithm could not find a solution
    cond = np.isnan(ubv_to_rubin_grid_filled)
    ubv_to_rubin_grid_final = np.zeros_like(ubv_to_rubin_grid_filled)
    ubv_to_rubin_grid_final[cond] = ubv_to_rubin_grid_nearest[cond]
    ubv_to_rubin_grid_final[~cond] = ubv_to_rubin_grid_filled[~cond]

    # Save the data
    grid_arr = np.squeeze(np.dstack([xx, yy]), axis=0)
    data_dir = '%s/data' % os.path.dirname(inspect.getfile(generate_ubv_to_rubin_grid))
    print(data_dir)
    ubv_to_rubin_filename = '%s/ubv_to_rubin-%s_grid.npz' % (data_dir, filter_name)
    print(ubv_to_rubin_filename)
    np.savez(ubv_to_rubin_filename,
             ubv_to_rubin_grid=ubv_to_rubin_grid_final.astype(np.float32),
             kdtree_grid=grid_arr.astype(np.float32))
    
def load_ubv_to_rubin_grid(filter_name):
    # Check for correct filter
    if filter_name not in ['u','g', 'r', 'i', 'z', 'y']:
        raise Exception("filter_name must be in: ['u','g', 'r', 'i', 'z', 'y']")

    # Load the ubv_to_rubin_grid from the file
    data_dir = '%s/data' % os.path.dirname(inspect.getfile(load_ubv_to_rubin_grid))
    ubv_to_rubin_filename = '%s/ubv_to_rubin-%s_grid.npz' % (data_dir, filter_name)
    ubv_to_rubin_grid_file = np.load(ubv_to_rubin_filename)

    # Generate a kdtree at the locations of all of the grid points
    ubv_to_rubin_grid = ubv_to_rubin_grid_file['ubv_to_rubin_grid']
    kdtree = cKDTree(ubv_to_rubin_grid_file['kdtree_grid'])

    return ubv_to_rubin_grid, kdtree

def transform_ubv_to_rubin(filter_name, ubv_B, ubv_V, ubv_R, ubv_U=None, ubv_I=None, twomass_J=None):
    """
    Converts ubv filters into rubin filters.

    Function is valid for values of:
        For filter u:
            0 < ubv_B - ubv_V < 6
            0 < ubv_U - ubv_B < 6
        
        For filters g, r, i, and z:
            0 < ubv_V - ubv_R < 6
            0 < ubv_B - ubv_V < 6
            
        For filter y:
            0 < ubv_I - twomass_J < 6
            0 < ubv_R - ubv_I < 6
    Parameters
    ----------
    filter_name : str
        rubin filter name of converted photometry
        Can be 'g', 'r', or 'i'
        If converting to rubin_i, then ubv_I must be provided

    ubv_B : array of floats
        ubv_B photometry of galaxia / SPISEA sources

    ubv_V : array of floats
        ubv_V photometry of galaxia / SPISEA sources

    ubv_R : array of floats
        ubv_R photometry of galaxia / SPISEA sources

    ubv_I : array of floats, optional
        ubv_I photometry of galaxia / SPISEA sources,
        required for converting to rubin_i

    Returns
    -------
    rubin_mag : array of floats
        rubin photometry of galaxia / SPISEA sources in `filter_name`

    """
    # Check for correct filter
    if filter_name not in ['u','g', 'r', 'i', 'z', 'y']:
        raise Exception("filter_name must be in: ['u','g', 'r', 'i', 'z', 'y']")

    if filter_name == 'i' and ubv_I is None:
        raise Exception('ubv_I must be provided to convert to rubin_i')
    elif filter_name == 'z' and ubv_I is None:
        raise Exception('ubv_I must be provided to convert to rubin_z')
    elif filter_name == 'u' and ubv_U is None:
        raise Exception('ubv_U must be provided to convert to rubin_u')
    elif filter_name == 'y' and ubv_I is None:
        if twomass_J is None:
            raise Exception('ubv_I and twomass_J must be provided to convert to rubin_y')
        

    # Convert the ubv photometry into the right format
    if filter_name == 'u':
        x_data = ubv_B - ubv_V
        y_data = ubv_U - ubv_B
    elif filter_name == 'y':
        x_data = ubv_I - twomass_J
        y_data = ubv_R - ubv_I
    else:
        x_data = ubv_V - ubv_R
        y_data = ubv_B- ubv_V

    data = np.squeeze(np.dstack([x_data, y_data]), axis=0)

    # Only query on data that is luminous
    cond_lum = ~np.isnan(data).any(axis=1)

    # Start with an empty array of nans
    rubin_diff = np.full(len(ubv_B), np.nan)

    # Find locations on the grid where x_data and y_data are located.
    # Put those values into the rubin_diff arrays
    ubv_to_rubin_grid, kdtree = load_ubv_to_rubin_grid(filter_name)
    _, indexes = kdtree.query(data[cond_lum])
    rubin_diff[cond_lum] = ubv_to_rubin_grid.flatten()[indexes]

    # Convert to rubin_g, rubin_r or rubin_i
    if filter_name == 'u':
        rubin_mag = ubv_U + rubin_diff
    elif filter_name == 'g':
        rubin_mag = ubv_V + rubin_diff
    elif filter_name == 'r':
        rubin_mag = ubv_R + rubin_diff
    elif filter_name == 'i':
        rubin_mag = ubv_I + rubin_diff
    elif filter_name == 'z':
        rubin_mag = ubv_I + rubin_diff
    elif filter_name == 'y':
        rubin_mag = twomass_J + rubin_diff

    return rubin_mag