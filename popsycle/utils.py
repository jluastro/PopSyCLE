#! /usr/bin/env python
"""
utils.py
Functions (and their associated functions) for utilities.
"""
import subprocess
import numpy as np
from astropy.table import Table, MaskedColumn
from astropy.table import vstack
import os
import gc
from popsycle import synthetic, orbits
import copy
import time
import datetime


def add_precision64(input_array, power):
    """
    Need more precision for kdtree to run properly. Convert inputs from
    float32 to float64, and add a random perturbation beginning in the
    nths decimal place.

    Parameters
    ----------
    input_array : float or array (float32)
        Thing that needs more precision.

    power : float
        To what place you want the perturbation.

    Returns
    -------
    output_array : float or array (float64)
        Thing that has more precision.

    """
    input_array = np.atleast_1d(input_array)

    # Perturb.
    pert = 10 ** power * (np.random.rand(len(input_array)) - 0.5)

    # Convert to float64.
    output_array = np.atleast_1d(np.float64(input_array))

    # Add the perturbation.
    output_array = output_array + pert

    return output_array


def make_symlinks(input_root, output_root):
    """
    Makes symlinks for galaxia_params.txt and
    perform_pop_syn.log files which are required for refine_events()
    
    Parameters
    ----------
    input_root : str
        Input root of galaxia_params.txt and perform_pop_syn.log.
        (Assumes they have the same input root).
    
    output_root : str
        Output root of galaxia_params.txt and perform_pop_syn.log.
    """
    links = ['_galaxia_params.txt', '_perform_pop_syn.log']
    for link in links:
        try:
            os.symlink(input_root + link, output_root + link)
        except:
            continue

def sample_spherical(npoints, speed, ndim=3):
    """
    Randomly sample points on a sphere.
    I found this code on stackexchange.

    Parameters
    ----------
    npoints : float
        The number of points you want to generate.

    speed : float
        The radius of the sphere (aka the magnitude of the vectors.)

    dim : float
        The dimension of the space in which the sphere is embedded
        (ndim = 3 samples points on a 2-sphere, aka a "normal" sphere)

    Returns
    -------
    An array of the vectors.
    """
    # Check that the speed vector is either a float or an array of length npoints
    if type(speed) != int:
        if type(speed) != float:
            if len(speed) != npoints:
                raise ValueError("{speed} must be either an int, float "
                                 "or array of length {npoints}")

    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    vec *= speed
    return vec


def return_nearest_gridpoint(grid, x_grid_arr, y_grid_arr, x_data, y_data):
    """
    Algorithm for finding the nearest grid cell on a 2D array given a
    datapoint that falls within the bounds of the 2D array.

    Parameters
    ----------
    grid : 2D numpy array
        2D array with size (len(y_grid_arr), len(x_grid_array))

    x_grid_arr : numpy array
        2D grid indices in the x-dimension

    y_grid_arr : numpy array
        2D grid indices in the y-dimension

    x_data : numpy array
        x-coordinate for data that will be located onto the grid

    y_data : numpy array
        y-coordinate for data that will be located onto the grid

    Returns
    -------
    gridpoint_arr : numpy array
        list of nearest cell values on the grid at
        the location of (x_data, y_data)

    """
    # Convert x_data and y_data to array if single data point is received
    x_data = np.atleast_1d(x_data)
    y_data = np.atleast_1d(y_data)

    gridpoint_arr = []
    for x, y in zip(x_data, y_data):
        # Loop through x_data and y_data
        if np.isnan(x) or np.isnan(y):
            # If either x_data or y_data is nan, return nan
            gridpoint = np.nan
        else:
            # Find location on the grid where x_data and y_data are
            # closest to the grid indices
            x_idx = np.argmin(np.abs(x - x_grid_arr))
            y_idx = np.argmin(np.abs(y - y_grid_arr))
            gridpoint = grid[y_idx, x_idx]
        gridpoint_arr.append(gridpoint)

    # Convert gridpoint_arr into numpy array
    gridpoint_arr = np.array(gridpoint_arr)

    # If only a single data point was received, return a single value
    if len(gridpoint_arr) == 1:
        gridpoint_arr = gridpoint_arr[0]

    return gridpoint_arr


def execute(cmd, shell=False):
    """
    Executes a command line instruction, captures the stdout and stderr

    Parameters
    ----------
    cmd : str
        Command line instruction, including any executables and parameters

    shell : bool, optional
        Determines if the command is run through the shell. Default is False.

    Returns
    -------
    stdout : str
        Contains the standard output of the executed process

    stderr : str
        Contains the standard error of the executed process

    """
    # Split the argument into a list suitable for Popen
    args = cmd.split()
    # subprocess.PIPE indicates that a pipe
    # to the standard stream should be opened.
    process = subprocess.Popen(args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               shell=shell)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')

    return stdout, stderr


def calc_normalized_counts(mag):
    """
    Parameters
    ----------
    mag: float
    
    Returns
    -------
    counts: int
    
    CHANGE THIS CODE IN THE FUTURE TO TAKE IN DIFFERENT ZERO POINTS! NEEDS ADDITIONAL DOCUMENTATION
    Right now this is only applicable for OGLE I band
    See confluence wiki for where these values come from...
    """
    a = -0.67815949
    b = 15.37993393
    counts = 10 ** (b + a * mag)
    return counts

def combine_refined_events(refined_events_filenames, overwrite=False,
                           output_file='default'):
    """
    Creates a combined refined_events out of multiple refined_events files

    Parameters
    ----------
    refined_events_filenames : list of strs
        Filenames of refined_events tables to be combined

    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.
        Default is False.

    output_file : str, optional
        The name of the final refined_events file.
        If set to 'default', the format will be generated from the first
        filename in refined_events_filenames following:
        combined_refined_events_<photometric_system>_<filt>_<red_law>.fits
    """
    # Check to see that all refined_events tables are unique
    if len(refined_events_filenames) != len(set(refined_events_filenames)):
        raise Exception('Duplicate filename found in '
                        'refined_events_filenames. Exiting...')

    # Loop over filenames, checking that each one exists
    print('Creating combined refined_events')
    refined_events_arr = []
    for filename in refined_events_filenames:
        if not os.path.exists(filename):
            raise Exception(f'{filename} cannot be found. Skipping...')
        print(f'-- Loading {filename}')
        refined_events = Table.read(filename)
        refined_events['original_filename'] = filename
        refined_events_arr.append(refined_events)

    # Combine astropy tables
    combined_refined_events = vstack(refined_events_arr)

    # Create an output_filename following the first filename
    if output_file == 'default':
        base = refined_events_filenames[0]
        output_file = 'combined_'
        output_file += '_'.join(base.split('_')[-5:])

    # Overwrite any exiting file with the same name
    if overwrite and os.path.exists(output_file):
        os.remove(output_file)

    # Save combined table to disk
    print(f'Saving to {output_file}')
    Table(combined_refined_events).write(output_file)

    return output_file

def calc_delta_c_LL(fratio, u0, thetaE):
    """
    Calculate the maximum-ish centroid shift for a luminous
    lens, no neighbors

    Parameters
    ----------
    fratio : flux ratio of the lens to source, i.e. f_L/f_S

    u0 : impact parameter

    thetaE : Einstein radius

    """
    final_delta_array = np.zeros(len(u0))
    final_u_array = np.zeros(len(u0))
    for j in np.arange(len(u0)):
        f = fratio[j]
        u_array = np.linspace(0, u0[j], 20)
        delta_array = np.zeros(len(u_array))
        for i in np.arange(len(u_array)):
            u = u_array[i]
            delta_array[i] = (u - f * u ** 2 * np.sqrt(u ** 2 + 4)) / (2 + u ** 2 + f * u * np.sqrt(u ** 2 + 4)) + (u * f) / (1 + f)
        max_idx = np.argmax(delta_array)
        final_delta_array[j] = delta_array[max_idx]
        final_u_array[j] = u_array[max_idx]
    final_delta_array = final_delta_array * thetaE

    return final_delta_array, final_u_array

def max_delta_c(u0, thetaE):
    max_delta_c_array = np.zeros(len(u0))

    big_idx = np.where(u0 > np.sqrt(2))[0]
    small_idx = np.where(u0 <= np.sqrt(2))[0]

    max_delta_c_array[big_idx] = calc_delta_c(u0[big_idx], thetaE[big_idx])
    max_delta_c_array[small_idx] = calc_delta_c(np.ones(len(small_idx)) * np.sqrt(2),
                                                thetaE[small_idx])

    return max_delta_c_array

def get_u_from_t(u0, t0, tE, t):
    """
    Given the time and separation at closest approach of lens and source
    and the Einstein radius, calculate the separation as a function of time.

    NOTE 1: You need to be consistent with your units for t0, tE, and t,
    i.e. pick whatever you want (days, years, etc.) but be self consistent.

    NOTE 2: There is a positive and negative solution for u.
    We return the positive solution.

    Parameters
    ----------
    u0 : float or array
        Minimum separation of lens and source (normalized to Einstein radius)

    t0 : float or array
        Time of minimum separation of lens and source

    tE : float or array
        Einstein crossing time of microlensing event

    t : float or array
        Time at which you want to calculate the separation

    Returns
    -------
    u : float or array
        Separation of lens and source (normalized to Einstein radius)
    """
    tau = (t - t0) / tE
    u = np.sqrt(u0 ** 2 + tau ** 2)
    return u


def get_t_from_u(u0, t0, tE, u):
    """
    Given the time and separation at closest approach of lens and source
    and the Einstein radius, calculate the time as a function of separation.

    NOTE 1: You need to be consistent with your units for t0, tE, and t,
    i.e. pick whatever you want (days, years, etc.) but be self consistent.

    NOTE 2: There is a positive and negative solution for t.
    We return the positive solution.

    Parameters
    ----------
    u0 : float or array
        Minimum separation of lens and source (normalized to Einstein radius)

    t0 : float or array
        Time of minimum separation of lens and source

    tE : float or array
        Einstein crossing time of microlensing event

    u : float or array
        Separation of lens and source (normalized to Einstein radius)

    Returns
    -------
    t : float or array
        Time corresponding to the separation u
    """
    t = tE * np.sqrt(u ** 2 - u0 ** 2) + t0
    return t


def calc_new_position(l0, b0, mu_lcosb, mu_b, t):
    """
    Given an initial position and proper motions in l and b,
    calculate the new position at some later time.

    Parameters
    ----------
    l0 : float or array
        Initial longitude, in DEGREES

    b0 : float or array
        Initial latitude, in DEGREES

    mu_lcosb : float or array
        Longitudinal proper motion l * cos(b), in MAS/YEAR

    mu_b : float or array
        Latitudinal roper motion, in MAS/YEAR

    t : float or array
        Time, in DAYS

    Returns
    -------
    l : float or array
        Latitude, in DEGREES

    b : float or array
        Longitude, in DEGREES
    """

    cos_b0 = np.cos(np.radians(b0))

    l = l0 + t * mu_lcosb * masyr_to_degday / cos_b0
    b = b0 + t * mu_b * masyr_to_degday

    return l, b


def calc_centroid_shift(glat_S, glon_S, glat_N, glon_N, f_L, f_S, f_N, u):
    """
    Calculate the centroid (astrometric) shift
    for a luminous lens and neighbors

    Parameters
    ----------
    glat_x : float or array
        Longitude of x (L = lens, S = source, N = neighbor centroid)

    glon_x : float or array
        Latitude of x (L = lens, S = source, N = neighbor centroid)

    f_x : float or array
        Flux of x (L = lens, S = source, N = all neighbors)

    u : float or array
        Dimensionless separation

    Returns
    -------
    delta_c_obs : float or array
        Magnitude of observed astrometric shift, in mas
    """
    ##########
    # Calculating the centroid shift in the frame of the lens
    ##########

    t1a1_t2a2 = (u ** 2 + 3) / (u * np.sqrt(u ** 2 + 4))
    a1_a2 = calc_magnification(u)

    glat_c_lensed = (t1a1_t2a2 * f_S * glat_S + glat_N * f_N) / (
            a1_a2 * f_S + f_L + f_N)
    glon_c_lensed = (t1a1_t2a2 * f_S * glon_S + glon_N * f_N) / (
            a1_a2 * f_S + f_L + f_N)

    glat_c_unlensed = (glat_S * f_S + glat_N * f_N) / (f_S + f_L + f_N)
    glon_c_unlensed = (glon_S * f_S + glon_N * f_N) / (f_S + f_L + f_N)

    glat_delta_c = glat_c_lensed - glat_c_unlensed
    glon_delta_c = glon_c_lensed - glon_c_unlensed

    delta_c_obs = np.sqrt(glat_delta_c ** 2 + glon_delta_c ** 2 * np.cos(
        np.radians(glat_S)) ** 2)

    # Convert from degrees to mas
    delta_c_obs *= 3.6 * 10 ** 6

    return delta_c_obs


def calc_f(lambda_eff):
    """
    Calculate that coefficient f that multiples E(B-V) to get the
    extinction in magnitudes
    """
    B = get_Alambda_AKs('Damineli16', 0.445)
    V = get_Alambda_AKs('Damineli16', 0.551)
    L = get_Alambda_AKs('Damineli16', lambda_eff)

    f = L * (B - V) ** -1

    return f


def generate_isochrones(iso_dir, logAge_min=5.01, logAge_max=10.291,
                        logAge_delta=0.01, metallicity=0.0,
                        include_ztf=True):
    """
    Generates isochrones needed for PopSyCLE simulations

    Parameters
    ----------
    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)

    logAge_min : float
        Minimum log10(age) of generated isochrones

    logAge_max : float
        Maximum log10(age) of generated isochrones

    logAge_delta : float
        Resolution in log10(age) of generated isochrones

    metallicity : float
        The metallicity in [M/H]

    include_ztf : bool, optional
        Determines if ztf filters should be included
        in the generated isochrones. Default is True.

    """
    from popsycle.synthetic import all_filt_list 
    from spisea import evolution, synthetic
    my_filt_list = all_filt_list
    # Add ztf filters to generated isochrones
    if include_ztf:
        my_filt_list += ['ztf,g', 'ztf,r', 'ztf,i']

    # Enforce sensible values
    logAge_min = max(5.01, logAge_min)
    logAge_max = min(10.291, logAge_max)
    logAge_delta = max(0.01, logAge_delta)

    # All possibles ages to two digit precision
    log_age_arr = np.arange(logAge_min, logAge_max, logAge_delta)
    print('Generating %i isochrones for '
          '%.2f < log_age < %.2f' % (len(log_age_arr),
                                     log_age_arr[0],
                                     log_age_arr[-1]))
    for i, log_age in enumerate(log_age_arr):
        print('Generating log_age = %.2f (%03d / %03d)' % (log_age, i+1, len(log_age_arr)))
        my_iso = synthetic.IsochronePhot(log_age, 0, 10,
                                         metallicity=metallicity,
                                         evo_model=evolution.MISTv1(),
                                         filters=my_filt_list,
                                         iso_dir=iso_dir)

        # Check that the isochrone has all of the filters in my_filt_list
        # If not, force recreating the isochrone with recomp=True
        my_iso_filters = [f for f in my_iso.points.colnames if 'm_' in f]
        my_filt_list_fmt = ['m_%s' % f.replace(',', '_') for f in my_filt_list]
        if set(my_filt_list_fmt) != set(my_iso_filters):
            _ = synthetic.IsochronePhot(log_age, 0, 10,
                                        metallicity=metallicity,
                                        evo_model=evolution.MISTv1(),
                                        filters=my_filt_list,
                                        iso_dir=iso_dir,
                                        recomp=True)
       
            
def remove_nan_companions(event_table_loc, comp_table_loc, event_output_loc = None, comp_output_loc = None, overwrite = False, rewrite_old_loc = None):
    """
    utility to remove companions that were not assigned a final mass
    or other parameters. This functionality will be implemented in
    perform_pop_syn, but this modifies those that were run with the
    nan companions included.
    
    If a SPISEA zams_mass of a companion is too small,
    the isochrone interpolator may not be able to assign a mass
    and other parameters so they'll be nan.
    This is especially a problem for old (>~ 9 Gyr) and high
    metallicity (Fe/H>0.4) objects.
    
    Parameters
    ----------
    event_table_loc : string
        Location of a companion table to modify
        
    comp_table_loc : string
        Location of a companion table to modify
    
    event_output_loc : string or None, optional
        String of output location of companion table without nan companions
        formated as '/path/to/file.fits'.
        If none, gets saved in the input location with '_no_nan' suffix
        Default is None.
        
    comp_output_loc : string or None, optional
        String of output location of companion table without nan companions
        formated as '/path/to/file.fits'.
        If none, gets saved in the input location with '_no_nan' suffix
        Default is None.
    
    overwrite : bool, optional
        Whether file at output_loc should be overwritten.
        Default is False.
        
    rewrite_old_loc : string or None, optional
        String of output location of companion table without nan companions.
        If None, does not get saved.
        Default is None.
    """
    t0 = time.time()
    now = datetime.datetime.now()
    
    event_table = Table.read(event_table_loc)
    comp_table = Table.read(comp_table_loc).filled(np.nan)
    original_comp_len = len(comp_table)
    if rewrite_old_loc is not None:
        rewrite_old_loc = comp_table.write(rewrite_old_loc)
    
    # Find companions with nan parameters
    bad_idxs = np.where(np.isnan(comp_table['mass']) == True)[0].tolist()
        
    # Set primary to not have companion
    if type(event_table['systemMass_L']) == MaskedColumn:
        event_table['systemMass_L'] = event_table['systemMass_L'].filled(np.nan)
    if type(event_table['systemMass_S']) == MaskedColumn:
        event_table['systemMass_S'] = event_table['systemMass_S'].filled(np.nan)
    for idx in bad_idxs:
        prim_type = comp_table['prim_type'][idx]
        system_idx = comp_table['system_idx'][idx]
        event_rows = event_table['obj_id_{}'.format(prim_type)] == system_idx
        # Decreases number of companions by 1 and if none are left, says system is no longer multiple
        event_table['N_companions_{}'.format(prim_type)][event_rows] -= 1
        companionless_event_rows = (event_table['N_companions_{}'.format(prim_type)] == 0) & (event_rows)
        event_table['isMultiple_{}'.format(prim_type)][companionless_event_rows] = 0

        # Fix system mass for those with dropped companions
        # Luminosity doesn't need to be fixed since nan mag counted as 0 flux
        # For those with no more companions, set systemMass to mass of primary
        systemMass_new_mass = event_table['mass_{}'.format(prim_type)][companionless_event_rows]
        event_table['systemMass_{}'.format(prim_type)][companionless_event_rows] = systemMass_new_mass 
        # For those with remaining companions, recalculate systemMass
        more_companions_event_rows = (event_table['N_companions_{}'.format(prim_type)] > 0) & (event_rows)
        systemMass_new =  event_table['mass_{}'.format(prim_type)][more_companions_event_rows]
        companions_mass = np.nan_to_num(comp_table['mass'][comp_table['system_idx'] == system_idx]) #treats nan masses as 0
        systemMass_new += sum(companions_mass) # adds the same to all of them, since it'd be multiple events with the same system if there are multiple rows
        event_table['systemMass_{}'.format(prim_type)][more_companions_event_rows] = systemMass_new

        # Fix Period, alpha, and phi
        # For those with remaining companions
        remaining_companions_rows = (comp_table['system_idx'] == system_idx) & (np.isnan(comp_table['mass']) == False)
        if sum(remaining_companions_rows) > 0:
            comp_table['P'][remaining_companions_rows] = orbits.a_to_P(systemMass_new, 10**comp_table['log_a'][remaining_companions_rows])
            event_table_df = event_table[more_companions_event_rows].to_pandas()
            event_table_df = event_table_df.set_index(['obj_id_L', 'obj_id_S'])
            comp_table_df = comp_table[remaining_companions_rows].to_pandas()
            comp_table_df = comp_table_df.set_index(['obj_id_L', 'obj_id_S'])
            joined_table = comp_table_df.join(event_table_df, lsuffix='_comp', rsuffix='_prim')
            alphas, phi_pi_Es, phis = synthetic.calculate_binary_angles(joined_table)
            comp_table['alpha'][remaining_companions_rows] = alphas
            comp_table['phi'][remaining_companions_rows] = phis
            
    # Remove bad companions
    comp_table.remove_rows(bad_idxs)

    if event_output_loc is None:
        event_output_loc = event_table_loc[0:-5] + '_no_nan.fits'
        
    if comp_output_loc is None:
        comp_output_loc = comp_table_loc[0:-5] + '_no_nan.fits'
    event_table.write(event_output_loc, overwrite=overwrite)
    comp_table.write(comp_output_loc, overwrite=overwrite)
    
    t1 = time.time()
    
    # Make logfile
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'INPUT FILES' + '\n'
    line1 = 'event table , ' + str(event_table_loc) + '\n'
    line2 = 'companion table , ' + str(comp_table_loc) + '\n'
    
    line3 = 'COMPANION INFORMATION' + '\n'
    line4 = '{} out of {} companions trimmed'.format(len(bad_idxs), original_comp_len)
    
    line5 = 'OUTPUT FILES' + '\n'
    line6 = 'new events table, ' + str(event_output_loc) + '\n'
    line6b = 'new companion table, ' + str(comp_output_loc) + '\n'
    line7 = ''
    if rewrite_old_loc is not None:
        line7 = 'rewrote old companion table, ' + str(rewrite_old_loc) + '\n'

    line8 = 'OTHER INFORMATION' + '\n'
    line9 = str(t1 - t0) + ' : total runtime (s)' + '\n'
    line10 = str(now) + ': time started' + '\n'

    with open(event_table_loc + '_remove_nan.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, dash_line, 
                        line4, empty_line, line5, dash_line, line6, line6b, line7, 
                        empty_line, line8, dash_line, line9, line10])
    
    return


def nan_to_zero(flux):
    """
    Modify the input flux array (or any numpy array for that matter) to replace all
    nan values with zeros.

    Parameters
    ----------
    flux : numpy array or Column
        Array of flux values. Can be me a regular numpy array or a Masked array or column.

    Returns
    -------
    flux : numpy array or Column
        copy of array now with nan's replaced with zeros
    """
    if type(flux) == np.ma.core.MaskedArray or type(flux) == MaskedColumn:
        flux = np.nan_to_num(flux.filled(np.nan))
    else:
        flux = np.nan_to_num(flux)

    return flux
