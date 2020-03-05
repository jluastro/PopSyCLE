#! /usr/bin/env python
"""
synthetic.py
Functions (and their associated functions) for running the PopSyCLE pipeline.
Including:
- write_galaxia_params
- perform_pop_syn
- calc_events
- refine_events
"""

import numpy as np
import h5py
import math
from astropy import units
from popsycle.filters import transform_ubv_to_ztf
from scipy.stats import maxwell
import astropy.coordinates as coord
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import Angle  # Angles
from astropy.table import Table
from astropy.table import vstack
from popstar.imf import imf
from popstar import synthetic, evolution, reddening, ifmr
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
import time
import datetime
from popsycle import ebf
import gc
import subprocess
import os
from sklearn import neighbors
import itertools
from multiprocessing import Pool
import inspect
import numpy.lib.recfunctions as rfn
from popsycle import utils


##########
# Conversions.
##########
masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcday = 1.0 * (3.086 * 10 ** 16) ** -1 * 86400.0
au_to_kpc = 4.848 * 10 ** -9

##########
# Prep converter function for determining
# the current-to-initial cluster mass vs. age.
# This global variable will get loaded the first time the function gets called.
##########
_Mclust_v_age_func = None

##########
# Dictionary for extinction law coefficients f_i, as a function of filter
# Damineli values from photometric bands (nm):
# B = 445, V = 551, I = 806, J = 1220, H = 1630, K = 2190, U = 365, R = 658
# Calculated using calc_f
# Schlegel and Schlafly photometric bands:
# B = 440, V = 543, I = 809, J = 1266, H = 1673, K = 2215, U = 337, R = 651
# ZTF photometric bands:
# G = 472.274, R = 633.961, I = 788.613
##########
filt_dict = {}
filt_dict['ubv_J'] = {'Schlafly11': 0.709, 'Schlegel99': 0.902, 'Damineli16': 0.662}
filt_dict['ubv_H'] = {'Schlafly11': 0.449, 'Schlegel99': 0.576, 'Damineli16': 0.344}
filt_dict['ubv_K'] = {'Schlafly11': 0.302, 'Schlegel99': 0.367, 'Damineli16': 0.172}
filt_dict['ubv_U'] = {'Schlafly11': 4.334, 'Schlegel99': 5.434, 'Damineli16': 5.022}
filt_dict['ubv_B'] = {'Schlafly11': 3.626, 'Schlegel99': 4.315, 'Damineli16': 3.757}
filt_dict['ubv_V'] = {'Schlafly11': 2.742, 'Schlegel99': 3.315, 'Damineli16': 2.757}
filt_dict['ubv_I'] = {'Schlafly11': 1.505, 'Schlegel99': 1.940, 'Damineli16': 1.496}
filt_dict['ubv_R'] = {'Schlafly11': 2.169, 'Schlegel99': 2.634, 'Damineli16': 2.102}
filt_dict['ztf_g'] = {'Damineli16': 3.453}
filt_dict['ztf_r'] = {'Damineli16': 2.228}
filt_dict['ztf_i'] = {'Damineli16': 1.553}

##########
# Dictionary for listing out supported photometric systems and filters
##########
photometric_system_dict = {}
photometric_system_dict['ubv'] = ['J', 'H', 'K', 'U', 'B', 'V', 'I', 'R']
photometric_system_dict['ztf'] = ['g', 'r']

##########
# List of all supported photometric systems and filters with PyPopStar labels
##########
all_filt_list = ['ubv,U', 'ubv,B', 'ubv,V', 'ubv,I', 'ubv,R',
                 'ukirt,H', 'ukirt,K', 'ukirt,J', 'ztf,g', 'ztf,r']

###########################################################################
############# Population synthesis and associated functions ###############
###########################################################################


def write_galaxia_params(output_root,
                         longitude, latitude, area,
                         seed=None):
    """
    Given an object root, sky location and area, creates the parameter
    file that Galaxia requires for running. User can also specify a seed for
    Galaxia to use in its object generation.

    Parameters
    ----------
    output_root : str
        The thing you want the output files to be named
        Examples:
           'myout'
           '/some/path/to/myout'
           '../back/to/some/path/myout'

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    Optional Parameters
    -------------------
    seed : int
         Seed Galaxia will use to generate objects. If not set, script will
         generate a seed from the current time. Setting this seed guarantees
         identical results.

    Outputs
    -------
    galaxia_params.<output_root>.txt : text file
        A text file with the parameters that Galaxia requires to run.
    """

    if seed is None:
        # Grab last two digits of the current Epochal time
        seed = int(str(time.time())[-2:])

    params = [
        "outputFile %s" % output_root,
        "outputDir ./",
        "photoSys UBV",
        "magcolorNames V,B-V",
        "appMagLimits[0] -1000",
        "appMagLimits[1] 1000",
        "absMagLimits[0] -1000",
        "absMagLimits[1] 1000",
        "colorLimits[0] -1000",
        "colorLimits[1] 1000",
        "geometryOption 1",
        "longitude %f" % longitude,
        "latitude %f" % latitude,
        "surveyArea %.5f" % area,
        "fSample 1",
        "popID -1",
        "warpFlareOn 1",
        "seed %i" % seed,
        "r_max 30",
        "starType 0",
        "photoError 0"
    ]

    galaxia_param_fname = 'galaxia_params.%s.txt' % output_root

    print('** Generating %s **' % galaxia_param_fname)

    with open(galaxia_param_fname, 'w') as f:
        for param in params:
            f.write(param + '\n')
            print('-- %s' % param)


def run_galaxia(output_root, longitude, latitude, area,
                seed=None):
    """
    Given an object root, sky location and area, creates the parameter
    file that Galaxia requires for running and executes Galaxia.
    User can also specify a seed for Galaxia to use in its object generation.

    Parameters
    ----------
    output_root : str
        The thing you want the output files to be named
        Examples:
           'myout'
           '/some/path/to/myout'
           '../back/to/some/path/myout'

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    Optional Parameters
    -------------------
    seed : int
         Seed Galaxia will use to generate objects. If not set, script will
         generate a seed from the current time. Setting this seed guarantees
         identical results.
    """

    # Writes out galaxia params to disk
    write_galaxia_params(output_root=output_root,
                         longitude=longitude,
                         latitude=latitude,
                         area=area,
                         seed=seed)

    # Execute Galaxia
    cmd = 'galaxia -r galaxia_params.%s.txt' % output_root
    print('** Executing Galaxia with galaxia_params.%s.txt **' % output_root)
    _ = utils.execute(cmd)
    print('** Galaxia complete **')


def perform_pop_syn(ebf_file, output_root, iso_dir,
                    bin_edges_number=None,
                    BH_kick_speed_mean=50, NS_kick_speed_mean=400,
                    additional_photometric_systems=None,
                    overwrite=False, seed=None):
    """
    Given some galaxia output, creates compact objects. Sorts the stars and
    compact objects into latitude/longitude bins, and saves them in an HDF5 file.

    Parameters
    ----------
    ebf_file : str or ebf file
        str : name of the ebf file from Galaxia
        ebf file : actually the ebf file from Galaxia

    output_root : str
        The thing you want the output files to be named
        Examples:
           'myout'
           '/some/path/to/myout'
           '../back/to/some/path/myout'

    iso_dir : filepath
        Where are the isochrones stored (for PopStar)

    Optional Parameters
    -------------------
    bin_edges_number : int
         Number of edges for the bins (bins = bin_edges_number - 1)
         Total number of bins is (bin_edges_number - 1)**2

    BH_kick_speed_mean : float
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.

    NS_kick_speed_mean : float
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract

    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / PyPopStar's ubv photometry and appended to the output files.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.
        Default is False.

    seed : int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for PyPopStar and PopSyCLE.
        Default None.

    Outputs
    -------
    <output_root>.h5 : hdf5 file
        NOTE: This is what _bin_lb_hdf5 returns.
        An hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects
        and stars sorted into those bins.

    <output_root>_label.fits : fits file
        NOTE: This is what make_label_file returns.
        A fits file that shows the correspondence between dataset name,
        latitude, longitude, and number of objects in that bin.
    """
    ##########
    # Error handling: check whether files exist and
    # whether input types are correct.
    ##########

    if not overwrite:
        # Check if HDF5 file exists already. If it does, throw an error message
        # to complain and exit.
        if os.path.isfile(output_root + '.h5'):
            raise Exception(
                'That .h5 file name is taken! Either delete the .h5 file, '
                'or pick a new name.')

        # Ditto with the fits label file.
        if os.path.isfile(output_root + '_label.fits'):
            raise Exception(
                'That .fits file name is taken! Either delete the .fits file, '
                'or pick a new name.')

    # Error handling/complaining if input types are not right.
    if ebf_file[-4:] != '.ebf':
        raise Exception('ebf_file must be an ebf file.')

    if type(output_root) != str:
        raise Exception('output_root must be a string.')

    if bin_edges_number is not None:
        if type(bin_edges_number) != int:
            raise Exception('bin_edges_number must be an integer.')

    if type(BH_kick_speed_mean) != int:
        if type(BH_kick_speed_mean) != float:
            raise Exception('BH_kick_speed_mean must be an integer or a float.')

    if type(NS_kick_speed_mean) != int:
        if type(NS_kick_speed_mean) != float:
            raise Exception('NS_kick_speed_mean must be an integer or a float.')

    if type(iso_dir) != str:
        raise Exception('iso_dir must be a string.')

    if seed is not None:
        if type(seed) != int:
            raise Exception('seed must be an integer.')

    if additional_photometric_systems is not None:
        if type(additional_photometric_systems) != list:
            raise Exception('additional_photometric_systems must either '
                            'None or a List (of strings).')

        for photometric_system in additional_photometric_systems:
            if photometric_system not in photometric_system_dict:
                raise Exception('strings in additional_photometric_systems '
                                'must be a valid option '
                                'in the photometric_system_dict.')

    ##########
    # Start of code
    #########

    # Set random seed
    np.random.seed(seed)

    t0 = time.time()

    #########
    # Read in only what you need
    # Control yourself... take only what you need from it
    # i.e. the log file, popid, and age
    ########
    t = ebf.read_ind(ebf_file, '/log', 0)
    popid_array = ebf.read(ebf_file, '/popid')
    age_array = ebf.read(ebf_file, '/age')

    # Just some counting/error checking things
    n_stars = len(popid_array)  # How many stars from Galaxia
    comp_counter = 0  # Number of compact objects made

    # Convert log to useful dictionary.
    ebf_log = make_ebf_log(t)

    ##########
    ###  Fetch the center point and area of the survey.
    ##########
    b = float(ebf_log['latitude'])
    l = float(ebf_log['longitude'])
    surveyArea = float(ebf_log['surveyArea'])

    ##########
    ### Make the bins for the latitude and longitude. All the information
    ### will be sorted into each of these bins in order to
    # handle large arrays, etc.
    ##########
    # Extend the edges a bit, that's what the * 1.1 is for
    # (to potentially catch any edge cases.)
    # make bins of size ~1/2 arcmin
    radius = np.sqrt(surveyArea / np.pi)  # degrees

    # Define bin_edges_number, if not given in input.
    if bin_edges_number is None:
        # st the widths are 1/2 arcmin
        bin_edges_number = int(60 * 2 * radius) + 1

    # Make sure we have enough bin edges (minimum is 3)
    if bin_edges_number < 2:
        bin_edges_number = 3
    lat_bin_edges = np.linspace(b - 1.1 * radius, b + 1.1 * radius,
                                bin_edges_number)
    long_bin_edges = np.linspace(l - 1.1 * radius, l + 1.1 * radius,
                                 bin_edges_number)
    # Angle wrapping for longitude
    wrap_id = np.where(long_bin_edges > 180)[0]
    long_bin_edges[wrap_id] -= 360

    ##########
    # Create h5py file to store lat/long binned output
    ##########
    h5file = h5py.File(output_root + '.h5', 'w')
    dataset = h5file.create_dataset('lat_bin_edges', data=lat_bin_edges)
    dataset = h5file.create_dataset('long_bin_edges', data=long_bin_edges)
    h5file.close()

    ##########
    # Reassign ages for stars that are less than logage 6
    # or greater than logage 10.01, since those are the limits of
    # PopStar evolution. Justification: in writeup/paper.
    ##########
    young_id = np.where(age_array <= 5.01)[0]
    age_array[young_id] = 5.0101
    old_id = np.where(age_array >= 10.14)[0]
    age_array[old_id] = 10.1399

    # Initialize a running counter for the "unique ID".
    next_id = n_stars  # for compact objects...
    n_binned_stars = 0  # for stars...

    ##########
    # Loop through population ID (i.e. bulge, disk, halo, ...) in order to
    # design bin sizes appropriate for each population. Dividing by population
    # is convenient because each has very different
    # age ranges and radius ranges.
    ##########
    for pid in range(10):
        print('*********************** Starting popid ' + str(pid))
        if np.sum(popid_array == pid) == 0:
            print('No stars with this pid. Skipping!')
            continue
        popid_idx = np.where(popid_array == pid)[0]
        if len(popid_idx) == 0:
            print('No stars with this pid. Skipping!')
            continue

        logage_min = np.min(age_array[popid_idx])
        logage_max = np.max(age_array[popid_idx])

        # For single age populations (popID = 7-9), have a single age bin...
        if logage_min == logage_max:
            logt_bins = np.array([logage_min * 0.99, logage_min * 1.01])

        # ...for multiple age populations (popID = 1-6,
        # break ages into 4 age bins.
        else:
            logt_bins = np.log10(np.logspace(logage_min, logage_max * 1.001, 5))

            # HARDCODED: Special handling for the youngest bin,
            # popID = 0 (Thin Disk < 150 Myr).
            if pid == 0:
                logt_bins = np.array([logage_min, 6.3, 7.0, 7.7, 8.0, 8.2])

        ##########
        # Loop through age bins
        #   -- note, don't slice tables as we don't want
        #   duplicated things carried in memory.
        ##########
        for aa in range(len(logt_bins) - 1):
            print('Starting age bin ', logt_bins[aa])
            # Mid-point age of bin.
            age_of_bin = (logt_bins[aa] + logt_bins[aa + 1]) / 2.0

            # Fetch the stars in this age bin.
            age_idx = np.where((age_array[popid_idx] >= logt_bins[aa]) &
                               (age_array[popid_idx] < logt_bins[aa + 1]))[0]
            len_adx = len(age_idx)

            # Figure out how many bins we will need.
            #   -- breaking up into managable chunks of 2 million stars each.
            num_stars_in_bin = 2e6
            num_bins = int(math.ceil(len_adx / num_stars_in_bin))

            ##########
            # Loop through bins of 2 million stars at a time.
            ##########
            for nn in range(num_bins):
                print('Starting sub-bin ', nn)
                n_start = int(nn * num_stars_in_bin)
                n_stop = int((nn + 1) * num_stars_in_bin)

                bin_idx = popid_idx[age_idx[n_start:n_stop]]

                ##########
                # Fill up star_dict
                ##########
                star_dict = {}
                star_dict['zams_mass'] = ebf.read_ind(ebf_file, '/smass', bin_idx)
                star_dict['mass'] = ebf.read_ind(ebf_file, '/mact', bin_idx)
                star_dict['px'] = ebf.read_ind(ebf_file, '/px', bin_idx)
                star_dict['py'] = ebf.read_ind(ebf_file, '/py', bin_idx)
                star_dict['pz'] = ebf.read_ind(ebf_file, '/pz', bin_idx)
                star_dict['vx'] = ebf.read_ind(ebf_file, '/vx', bin_idx)
                star_dict['vy'] = ebf.read_ind(ebf_file, '/vy', bin_idx)
                star_dict['vz'] = ebf.read_ind(ebf_file, '/vz', bin_idx)
                star_dict['age'] = age_array[bin_idx]
                star_dict['popid'] = popid_array[bin_idx]
                star_dict['exbv'] = ebf.read_ind(ebf_file, '/exbv_schlegel', bin_idx)
                star_dict['glat'] = ebf.read_ind(ebf_file, '/glat', bin_idx)
                star_dict['glon'] = ebf.read_ind(ebf_file, '/glon', bin_idx)
                star_dict['mbol'] = ebf.read_ind(ebf_file, '/lum', bin_idx)
                star_dict['grav'] = ebf.read_ind(ebf_file, '/grav', bin_idx)
                star_dict['teff'] = ebf.read_ind(ebf_file, '/teff', bin_idx)
                star_dict['feh'] = ebf.read_ind(ebf_file, '/feh', bin_idx)
                star_dict['rad'] = ebf.read_ind(ebf_file, '/rad', bin_idx)
                star_dict['rem_id'] = np.zeros(len(bin_idx))
                star_dict['obj_id'] = np.arange(len(bin_idx)) + n_binned_stars
                n_binned_stars += len(bin_idx)

                ##########
                # Angle wrapping for longitude
                ##########
                wrap_idx = np.where(star_dict['glon'] > 180)[0]
                star_dict['glon'][wrap_idx] -= 360

                ##########
                # Add UBV magnitudes
                ##########
                star_dict['ubv_J'] = ebf.read_ind(ebf_file, '/ubv_J', bin_idx)
                star_dict['ubv_H'] = ebf.read_ind(ebf_file, '/ubv_H', bin_idx)
                star_dict['ubv_K'] = ebf.read_ind(ebf_file, '/ubv_K', bin_idx)
                star_dict['ubv_U'] = ebf.read_ind(ebf_file, '/ubv_U', bin_idx)
                star_dict['ubv_I'] = ebf.read_ind(ebf_file, '/ubv_I', bin_idx)
                star_dict['ubv_B'] = ebf.read_ind(ebf_file, '/ubv_B', bin_idx)
                star_dict['ubv_V'] = ebf.read_ind(ebf_file, '/ubv_V', bin_idx)
                star_dict['ubv_R'] = ebf.read_ind(ebf_file, '/ubv_R', bin_idx)

                ##########
                # Add ztf magnitudes
                ##########
                if additional_photometric_systems is not None:
                    if 'ztf' in additional_photometric_systems:
                        # Pull out ubv magnitudes needed for photometric conversions
                        ubv_b = star_dict['ubv_B']
                        ubv_v = star_dict['ubv_V']
                        ubv_r = star_dict['ubv_R']
                        ubv_i = star_dict['ubv_I']

                    ztf_g = transform_ubv_to_ztf('g', ubv_b, ubv_v, ubv_r, ubv_i)
                    ztf_r = transform_ubv_to_ztf('r', ubv_b, ubv_v, ubv_r, ubv_i)
                    ztf_i = transform_ubv_to_ztf('i', ubv_b, ubv_v, ubv_r, ubv_i)
                    star_dict['ztf_g'] = ztf_g
                    star_dict['ztf_r'] = ztf_r
                    star_dict['ztf_i'] = ztf_i

                ##########
                # Add spherical velocities vr, mu_b, mu_lcosb
                ##########
                vr, mu_b, mu_lcosb = calc_sph_motion(star_dict['vx'],
                                                     star_dict['vy'],
                                                     star_dict['vz'],
                                                     star_dict['rad'],
                                                     star_dict['glat'],
                                                     star_dict['glon'])
                #########
                # Add precision to r, b, l, vr, mu_b, mu_lcosb
                #########
                star_dict['rad'] = utils.add_precision64(star_dict['rad'], -4)
                star_dict['glat'] = utils.add_precision64(star_dict['glat'], -4)
                star_dict['glon'] = utils.add_precision64(star_dict['glon'], -4)
                star_dict['vr'] = utils.add_precision64(vr, -4)
                star_dict['mu_b'] = utils.add_precision64(mu_b, -4)
                star_dict['mu_lcosb'] = utils.add_precision64(mu_lcosb, -4)

                ##########
                # Perform population synthesis.
                ##########
                mass_in_bin = np.sum(star_dict['zams_mass'])

                stars_in_bin = {}
                for key, val in star_dict.items():
                    stars_in_bin[key] = val

                comp_dict, next_id = _make_comp_dict(iso_dir, age_of_bin,
                                                     mass_in_bin, stars_in_bin, next_id,
                                                     BH_kick_speed_mean=BH_kick_speed_mean,
                                                     NS_kick_speed_mean=NS_kick_speed_mean,
                                                     additional_photometric_systems=additional_photometric_systems,
                                                     seed=seed)

                ##########
                #  Bin in l, b all stars and compact objects.
                ##########
                if comp_dict is not None:
                    comp_counter += len(comp_dict['mass'])
                    _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                                 comp_dict, output_root)
                _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                             stars_in_bin,
                             output_root)
                ##########
                # Done with galaxia output in dictionary t and ebf_log.
                # Garbage collect in order to save space.
                ##########
                del star_dict
                gc.collect()

    t1 = time.time()
    print('Total run time is {0:f} s'.format(t1 - t0))

    ##########
    # Figure out how much stuff got binned.
    ##########
    binned_counter = 0
    hf = h5py.File(output_root + '.h5', 'r')
    for key in hf:
        if 'bin_edges' not in key:
            binned_counter += len(hf[key])

    ##########
    # Make label file containing information about l,b bins
    ##########
    make_label_file(output_root, overwrite=overwrite)

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    microlens_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    microlens_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=microlens_path).decode('ascii').strip()
    popstar_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                           cwd=popstar_path).decode('ascii').strip()
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'ebf_file , ' + ebf_file + '\n'
    line2 = 'output_root , ' + output_root + '\n'
    line3 = 'bin_edges_number , ' + str(bin_edges_number) + '\n'
    line4 = 'BH_kick_speed_mean , ' + str(BH_kick_speed_mean) + ' , (km/s)' + '\n'
    line5 = 'NS_kick_speed_mean , ' + str(NS_kick_speed_mean) + ' , (km/s)' + '\n'
    line6 = 'iso_dir , ' + iso_dir + '\n'

    line7 = 'VERSION INFORMATION' + '\n'
    line8 = str(now) + ' : creation date' + '\n'
    line9 = popstar_hash + ' : PopStar commit' + '\n'
    line10 = microlens_hash + ' : microlens commit' + '\n'

    line11 = 'OTHER INFORMATION' + '\n'
    line12 = str(t1 - t0) + ' : total runtime (s)' + '\n'
    line13 = str(n_stars) + ' : total stars from Galaxia' + '\n'
    line14 = str(comp_counter) + ' : total compact objects made' + '\n'
    line15 = str(binned_counter) + ' : total things binned' + '\n'

    line16 = 'FILES CREATED' + '\n'
    line17 = output_root + '.h5 : HDF5 file' + '\n'
    line18 = output_root + '_label.fits : label file' + '\n'

    with open(output_root + '_perform_pop_syn.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, line4, line5,
                        line6, empty_line, line7, dash_line, line8, line9,
                        line10, empty_line, line11, dash_line, line12, line13,
                        line14, line15, empty_line, line16, dash_line, line17,
                        line18])

    ##########
    # Informative print statements.
    ##########
    if (n_stars + comp_counter) != binned_counter:
        print('***************** WARNING ******************')
        print('Number of things in != number of things out.')
        print('********************************************')

    print('******************** INFO **********************')
    print('Total number of stars from Galaxia: ' + str(n_stars))
    print('Total number of compact objects made: ' + str(comp_counter))
    print('Total number of things binned: ' + str(binned_counter))

    return


def calc_current_initial_ratio(iso_dir,
                               out_file='current_initial_stellar_mass_ratio.txt',
                               seed=None):
    """
    Makes 10**7 M_sun clusters in PopStar at various ages, to calculate the
    ratio of current to initial cluster mass. The range of ages goes
    from 6 to 10.089 log(age/yr). Creates a text file, the first column is
    the ages, and second column is the ratio.

    Parameters
    -------------------
    out_file : string
        Full name of the file to store the output columns:
        log(age/yr), current_to_initial_mass_ratio

    iso_dir : filepath
        Where you are storing isochrones

    Optional Parameters
    -------------------
    seed : int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for PyPopStar and PopSyCLE.
        Default None.

    Output
    ------
    out_file : txt file
         File of the initial-final cluster mass.

    """
    # Set up the input information for PopStar to make a cluster.
    # age in log(age/yr)
    logage_vec = np.concatenate((np.arange(5.01, 10.01 + 0.005, 0.2),
                                 np.array([10.04, 10.14])))
    cluster_mass = 10 ** 7  # mass in m_sol
    # IMF power law boundaries. changed start from 0.08 to 0.1 because of MIST.
    mass_limits = np.array([0.1, 0.5, 120])
    # IMF powers
    powers = np.array([-1.3, -2.3])
    # initialize current-initial ratio vector
    current_initial_ratio_array = np.zeros(len(logage_vec))

    # Run all the clusters for the different ages in logage_vec
    for i in range(len(logage_vec)):
        # make isochrone
        my_iso = synthetic.IsochronePhot(logage_vec[i], 0, 10,
                                         evo_model=evolution.MISTv1(),
                                         filters=all_filt_list,
                                         iso_dir=iso_dir)

        # define IMF
        trunc_kroupa = imf.IMF_broken_powerlaw(mass_limits, powers)
        # make cluster
        cluster = synthetic.ResolvedCluster(my_iso, trunc_kroupa, cluster_mass,
                                            seed=seed)
        output = cluster.star_systems

        # Find the stars in MIST not in Galaxia (i.e. WDs)
        # and figure out how much mass they contribute
        bad_idx = np.where(output['phase'] == 101)[0]
        bad_mass = np.sum(output['mass'][bad_idx])  # in m_sol

        # Figure out the current mass (stars + post-AGB stars)
        current_mass = np.sum(output['mass'])  # in m_sol

        # Figure out the current STELLAR mass (no post-AGB stars)
        # This is what we want for the ratio
        # since Galaxia does NOT have white dwarfs.
        current_stellar_mass = current_mass - bad_mass

        # Calculate the ratio of initial to final stellar mass.
        current_initial_mass_ratio = current_stellar_mass / cluster_mass
        current_initial_ratio_array[i] = current_initial_mass_ratio
        print(current_initial_mass_ratio)

    # Reshape vectors so they can be concatenated into two columns
    logage_vec = logage_vec.reshape((len(logage_vec), 1))
    current_initial_ratio_array = current_initial_ratio_array.reshape((len(logage_vec), 1))
    tab = np.concatenate((logage_vec, current_initial_ratio_array), axis=1)
    np.savetxt(out_file, tab)

    return


def current_initial_ratio(logage, ratio_file, iso_dir, seed=None):
    """
    Calculates the ratio of the current cluster mass to the initial
    mass of the cluster.

    Parameters
    ----------
    logage : float
        The age of the cluster, given in log(age/yr)

    ratio_file : txt file
        The name of the text file that contains the
        current-initial cluster mass ratio.

    iso_dir : filepath
        Where are the isochrones stored (for PopStar)

    Optional Parameters
    -------------------
    seed : int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for PyPopStar and PopSyCLE.
        Default None.

    Return
    ------
        The current-initial cluster mass ratio (float)

    """
    global _Mclust_v_age_func

    if _Mclust_v_age_func == None:
        try:
            boop = np.loadtxt(ratio_file)
        except Exception as e:
            calc_current_initial_ratio(iso_dir=iso_dir, out_file=ratio_file,
                                       seed=seed)
            boop = np.loadtxt(ratio_file)

        logage_vec = boop[:, 0]
        current_initial_mass_ratio = boop[:, 1]
        _Mclust_v_age_func = interp1d(logage_vec, current_initial_mass_ratio)

    return _Mclust_v_age_func(logage)


def _make_comp_dict(iso_dir, log_age, currentClusterMass, star_dict, next_id,
                    BH_kick_speed_mean=50, NS_kick_speed_mean=400,
                    additional_photometric_systems=None,
                    seed=None):
    """
    Perform population synthesis.

    Parameters
    ----------
    iso_dir : filepath
        Where are the isochrones stored (for PopStar)

    log_age : float
        log(age/yr) of the cluster you want to make

    currentClusterMass : float
        Mass of the cluster you want to make (M_sun)

    star_dict : dictionary
        The number of entries for each key is the number of stars.

    next_id : The next unique ID number (int) that will be assigned to
              the new compact objects created.

    Optional Parameters
    -------------------
    BH_kick_speed_mean : float
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.

    NS_kick_speed_mean : float
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract

    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / PyPopStar's ubv photometry and appended to the output files.

    seed : int
         Seed used to sample the kde tree. If set to any number,
         PyPopStar will also be forced to use 42 as a
         random seed for calls to ResolvedCluster.
         Default is None.

    Returns
    -------
    comp_dict : dictionary
        Keys are the same as star_dict, just for compact objects.

    next_id : int
        Updated next unique ID number (int) that will be assigned to
        the new compact objects created.

    """
    comp_dict = None

    # Calculate the initial cluster mass
    # changed from 0.08 to 0.1 at start because MIST can't handle.
    massLimits = np.array([0.1, 0.5, 120])
    powers = np.array([-1.3, -2.3])
    my_ifmr = ifmr.IFMR()
    ratio_file = '%s/current_initial_stellar_mass_ratio.txt' % iso_dir
    ratio = current_initial_ratio(logage=log_age,
                                  ratio_file=ratio_file,
                                  iso_dir=iso_dir,
                                  seed=seed)
    initialClusterMass = currentClusterMass / ratio

    ##########
    # Create the PopStar table (stars and compact objects).
    #    - it is only sensible to do this for a decent sized cluster.
    ##########
    if initialClusterMass > 100:
        # MAKE isochrone
        # -- arbitrarily chose AKs = 0, distance = 10 pc
        # (irrelevant, photometry not used)
        # Using MIST models to get white dwarfs
        my_iso = synthetic.IsochronePhot(log_age, 0, 10,
                                         evo_model=evolution.MISTv1(),
                                         filters=all_filt_list,
                                         iso_dir=iso_dir)

        # Check that the isochrone has all of the filters in filt_list
        # If not, force recreating the isochrone with recomp=True
        my_iso_filters = [f for f in my_iso.points.colnames if 'm_' in f]
        filt_list = ['m_%s' % f.replace(',', '_') for f in all_filt_list]
        if set(filt_list) != set(my_iso_filters):
            my_iso = synthetic.IsochronePhot(log_age, 0, 10,
                                             evo_model=evolution.MISTv1(),
                                             filters=all_filt_list,
                                             iso_dir=iso_dir,
                                             recomp=True)

        # !!! Keep trunc_kroupa out here !!! Death and destruction otherwise.
        # DON'T MOVE IT OUT!
        trunc_kroupa = imf.IMF_broken_powerlaw(massLimits, powers)

        # MAKE cluster
        cluster = synthetic.ResolvedCluster(my_iso, trunc_kroupa,
                                            initialClusterMass, ifmr=my_ifmr,
                                            seed=seed)
        output = cluster.star_systems

        # Create the PopStar table with just compact objects

        # The compact IDs are:
        # 101: WD, 102: NS, 103: BH
        # -1: PMS, 0: MS, 2: RGB, 3: CHeB, 4: EAGB, 5: TPAGB, 9: WR
        # WE CONSIDER EVERYTHING FROM THE MIST MODELS
        # WITH PHASE 6 (postAGB) TO BE WHITE DWARFS
        compact_ID = np.where((output['phase'] == 101) |
                              (output['phase'] == 102) |
                              (output['phase'] == 103))[0]
        comp_table = output[compact_ID]

        # Removes unused columns to conserve memory.
        keep_columns = ['mass', 'phase', 'mass_current', 'm_ubv_I', 'm_ubv_R',
                        'm_ubv_B', 'm_ubv_U', 'm_ubv_V', 'm_ukirt_H',
                        'm_ukirt_J', 'm_ukirt_K']
        if additional_photometric_systems is not None:
            if 'ztf' in additional_photometric_systems:
                keep_columns += ['m_ztf_g', 'm_ztf_r', 'm_ztf_i']
        comp_table.keep_columns(keep_columns)

        # Fill out the rest of comp_dict
        if len(comp_table['mass']) > 0:

            # Turn astropy table into dictionary to conserve memory.
            comp_dict = {}
            comp_dict['mass'] = comp_table['mass_current'].data
            comp_dict['rem_id'] = comp_table['phase'].data
            comp_dict['zams_mass'] = comp_table['mass'].data

            ##########
            # Assign spherical positions and velocities to all compact objects.
            ##########
            # Create velocities and positions using Kernel Density Estimation
            # Galactic position (r, l, b)
            # and heliocentric velocity (vx, vy, vz)
            kde_in_data = np.array([star_dict['rad'], star_dict['glat'],
                                    star_dict['glon'], star_dict['vx'],
                                    star_dict['vy'], star_dict['vz']]).T
            # .T because kde.fit needs the rows/columns this way.
            kde = neighbors.KernelDensity(bandwidth=0.0001)
            kde.fit(kde_in_data)

            kde_out_data = kde.sample(len(comp_dict['mass']),
                                      random_state=seed)
            comp_dict['rad'] = kde_out_data[:, 0]
            comp_dict['glat'] = kde_out_data[:, 1]
            comp_dict['glon'] = kde_out_data[:, 2]
            comp_dict['vx'] = kde_out_data[:, 3]
            comp_dict['vy'] = kde_out_data[:, 4]
            comp_dict['vz'] = kde_out_data[:, 5]

            ##########
            # Assign x, y, z position.
            ##########
            comp_helio = galactic_to_heliocentric(comp_dict['rad'],
                                                  comp_dict['glat'],
                                                  comp_dict['glon'])
            comp_dict['px'], comp_dict['py'], comp_dict['pz'] = comp_helio

            ##########
            # Add kicks to NSs and BHs.
            ##########

            # Maxwellian pdf is sqrt(2/pi) * x^2 * exp[-x^2/(2 * a^2)] / a^3,
            # with a scale parameter `a`.
            # The Maxwellian mean is 2 * a * sqrt(2/pi).
            # Here we calculate the scipy.stats `scale` by dividing the
            # user defined mean by the Maxwellian mean.

            NS_idx = np.where(comp_dict['rem_id'] == 102)[0]
            NS_kick_speed_scale = NS_kick_speed_mean / (2*np.sqrt(2/np.pi))
            if len(NS_idx) > 0:
                NS_kick_speed = maxwell.rvs(loc=0, scale=NS_kick_speed_scale, size=len(NS_idx))
                NS_kick = utils.sample_spherical(len(NS_idx), NS_kick_speed)
                comp_dict['vx'][NS_idx] += NS_kick[0]
                comp_dict['vy'][NS_idx] += NS_kick[1]
                comp_dict['vz'][NS_idx] += NS_kick[2]

            BH_idx = np.where(comp_dict['rem_id'] == 103)[0]
            BH_kick_speed_scale = BH_kick_speed_mean / (2*np.sqrt(2/np.pi))
            if len(BH_idx) > 0:
                BH_kick_speed = maxwell.rvs(loc=0, scale=BH_kick_speed_scale, size=len(BH_idx))
                BH_kick = utils.sample_spherical(len(BH_idx), BH_kick_speed)
                comp_dict['vx'][BH_idx] += BH_kick[0]
                comp_dict['vy'][BH_idx] += BH_kick[1]
                comp_dict['vz'][BH_idx] += BH_kick[2]

            # Add precision to r, b, l
            comp_dict['rad'] = utils.add_precision64(comp_dict['rad'], -4)
            comp_dict['glat'] = utils.add_precision64(comp_dict['glat'], -4)
            comp_dict['glon'] = utils.add_precision64(comp_dict['glon'], -4)

            # Assign vr, mu_b, mu_lcosb.
            comp_sph = calc_sph_motion(comp_dict['vx'],
                                       comp_dict['vy'],
                                       comp_dict['vz'],
                                       comp_dict['rad'],
                                       comp_dict['glat'],
                                       comp_dict['glon'])
            comp_dict['vr'], comp_dict['mu_b'], comp_dict['mu_lcosb'] = comp_sph

            # Add precision to vr, mu_b, mu_lcosb
            comp_dict['vr'] = utils.add_precision64(comp_dict['vr'], -4)
            comp_dict['mu_b'] = utils.add_precision64(comp_dict['mu_b'], -4)
            comp_dict['mu_lcosb'] = utils.add_precision64(comp_dict['mu_lcosb'], -4)

            # Assign age.
            comp_dict['age'] = log_age * np.ones(len(comp_dict['vx']))

            ##########
            # Initialize values for compact object photometry.
            # All of the values are np.nan
            # These are all the outputs from the IFMR of Raithel and Kalirai.
            ##########
            comp_dict['exbv'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_I'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_K'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_J'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_U'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_R'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_B'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_V'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['ubv_H'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['teff'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['grav'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['mbol'] = np.full(len(comp_dict['vx']), np.nan)
            comp_dict['feh'] = np.full(len(comp_dict['vx']), np.nan)

            if additional_photometric_systems is not None:
                if 'ztf' in additional_photometric_systems:
                    comp_dict['ztf_g'] = np.full(len(comp_dict['vx']), np.nan)
                    comp_dict['ztf_r'] = np.full(len(comp_dict['vx']), np.nan)
                    comp_dict['ztf_i'] = np.full(len(comp_dict['vx']), np.nan)

            ##########
            # FIX THE BAD PHOTOMETRY FOR LUMINOUS WHITE DWARFS
            # For non-dark WDs only (the ones from MIST):
            # Approximate extinction from the nearest (in 3-D space) star.
            # Get WD photometry from MIST models.
            ##########
            lum_WD_idx = np.argwhere(~np.isnan(comp_table['m_ubv_I']))

            if len(lum_WD_idx) > 0:
                star_xyz = np.array([star_dict['px'],
                                     star_dict['py'],
                                     star_dict['pz']]).T
                comp_xyz = np.array([comp_dict['px'][lum_WD_idx],
                                     comp_dict['py'][lum_WD_idx],
                                     comp_dict['pz'][lum_WD_idx]]).T

                kdt = cKDTree(star_xyz)
                dist, indices = kdt.query(comp_xyz)

                comp_dict['exbv'][lum_WD_idx] = star_dict['exbv'][indices.T]

                comp_dict['ubv_I'][lum_WD_idx] = comp_table['m_ubv_I'][lum_WD_idx].data
                comp_dict['ubv_K'][lum_WD_idx] = comp_table['m_ukirt_K'][lum_WD_idx].data
                comp_dict['ubv_J'][lum_WD_idx] = comp_table['m_ukirt_J'][lum_WD_idx].data
                comp_dict['ubv_U'][lum_WD_idx] = comp_table['m_ubv_U'][lum_WD_idx].data
                comp_dict['ubv_R'][lum_WD_idx] = comp_table['m_ubv_R'][lum_WD_idx].data
                comp_dict['ubv_B'][lum_WD_idx] = comp_table['m_ubv_B'][lum_WD_idx].data
                comp_dict['ubv_V'][lum_WD_idx] = comp_table['m_ubv_V'][lum_WD_idx].data
                comp_dict['ubv_H'][lum_WD_idx] = comp_table['m_ukirt_H'][lum_WD_idx].data
                if additional_photometric_systems is not None:
                    if 'ztf' in additional_photometric_systems:
                        comp_dict['ztf_g'][lum_WD_idx] = comp_table['m_ztf_g'][lum_WD_idx].data
                        comp_dict['ztf_r'][lum_WD_idx] = comp_table['m_ztf_r'][lum_WD_idx].data
                        comp_dict['ztf_i'][lum_WD_idx] = comp_table['m_ztf_i'][lum_WD_idx].data

                # Memory cleaning
                del comp_table
                gc.collect()

            # Assign population and object ID.
            comp_dict['popid'] = star_dict['popid'][0] * np.ones(len(comp_dict['vx']))
            comp_dict['obj_id'] = np.arange(len(comp_dict['vx'])) + next_id

            next_id += len(comp_dict['vx'])

    else:
        comp_dict = None

    return comp_dict, next_id


def _generate_comp_dtype(obj_arr):
    """
    Create compound datatype by looping over the keys of the obj_arr.
    Assigns integers as datatypes where reasonable, and float64 to the rest
    of the columns

    Parameters
    ----------
    obj_arr : array or None
        Array of stars or compact objects to be binned
        (either comp_dict or star_dict)

    Returns
    ------
    comp_dtype : np.dtype
        Numpy datatype with all of the keys of obj_arr
    """

    comp_dtype_arr = []
    for key in obj_arr:
        if key in ['rem_id', 'popid']:  # int16 (up to 32767)
            d = (key, 'i2')
        elif key in ['obj_id']:  # int32 (up to 2147483647)
            d = (key, 'i4')
        else:
            d = (key, 'f8')  # float64
        comp_dtype_arr.append(d)
    comp_dtype = np.dtype(comp_dtype_arr)
    return comp_dtype


def _bin_lb_hdf5(lat_bin_edges, long_bin_edges, obj_arr, output_root):
    """
    Given stars and compact objects, sort them into latitude and
    longitude bins. Save each latitude and longitude bin, and the edges that
    define the latitude and longitude bins, as datasets in a hdf5 file.

    Parameters
    ----------
    lat_bin_edges : array
        Edges for the latitude binning (deg)

    long_bin_edges : array
        Edges for the longitude binning (deg)

    obj_arr : array or None
        Array of stars or compact objects to be binned
        (either comp_dict or star_dict)

    output_root : str
        The path and name of the hdf5 file,
        without suffix (will be saved as output_root.h5)

    Optional Parameters
    -------------------

    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / PyPopStar's ubv photometry and appended to the output files.

    Output
    ------
    output_root.h5 : hdf5 file
        An hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects and stars sorted into
        those bins.
    """
    # Create compound datatype from obj_arr
    comp_dtype = _generate_comp_dtype(obj_arr)

    ##########
    # Loop through the latitude and longitude bins.
    ##########
    for ll in range(len(long_bin_edges) - 1):
        for bb in range(len(lat_bin_edges) - 1):
            # Open our HDF5 file for reading and appending.
            # Create as necessary.
            hf = h5py.File(output_root + '.h5', 'r+')

            # HDF5 dataset name
            dset_name = 'l' + str(ll) + 'b' + str(bb)

            # Create data set if needed. Start with 0 stars in the dataset.
            if dset_name not in hf:
                dataset = hf.create_dataset(dset_name, shape=(0,),
                                            chunks=(1e4,),
                                            maxshape=(None,),
                                            dtype=comp_dtype)
            else:
                dataset = hf[dset_name]

            ##########
            # Binning the stars and/or compact objects
            ##########
            if obj_arr is not None:
                id_lb = np.where((obj_arr['glat'] >= lat_bin_edges[bb]) &
                                 (obj_arr['glat'] < lat_bin_edges[bb + 1]) &
                                 (obj_arr['glon'] >= long_bin_edges[ll]) &
                                 (obj_arr['glon'] < long_bin_edges[ll + 1]))[0]

                if len(id_lb) == 0:
                    continue

                # Loop over the obj_arr and add all columns
                # (matching id_lb) into save_data
                save_data = np.empty(len(id_lb), dtype=comp_dtype)
                for colname in obj_arr:
                    save_data[colname] = obj_arr[colname][id_lb]

                # Resize the dataset and add data.
                old_size = dataset.shape[0]
                new_size = old_size + len(id_lb)
                dataset.resize((new_size, ))
                dataset[old_size:new_size] = save_data

            hf.close()

    return


############################################################################
########### Candidate event calculation and associated functions ###########
############################################################################

def calc_events(hdf5_file, output_root2,
                radius_cut=2, obs_time=1000, n_obs=101, theta_frac=2,
                blend_rad=0.65, n_proc=1,
                seed=None, overwrite=False):
    """
    Calculate microlensing events

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file.

    radius_cut : float
        Initial radius cut, in ARCSECONDS.

    obs_time : float
        Survey duration, in DAYS.

    n_obs : float
        Number of observations.

    theta_frac : float
        Another cut, in multiples of Einstein radii.

    blend_rad : float
        Stars within this distance of the lens are said to be blended.
        Units are in ARCSECONDS.

    Optional Parameters
    -------------------
    n_proc : int
        Number of processors to use. Should not exceed the number of cores.
        Default is one processor (no parallelization).

    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / PyPopStar's ubv photometry and appended to the output files.

    seed : int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for PyPopStar and PopSyCLE.
        Default None.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.


    Output
    ------
    <output_root2>_events.fits : Astropy .fits table
        Table of candidate microlensing events. The number of rows
        corresponds to the number of candidate events.

    """

    ##########
    # Error handling: check whether files exist and
    # whether input types are correct.
    ##########

    # Check if .fits file exists already. If it does, throw an error message
    # to complain and exit.
    if not overwrite:
        if os.path.isfile(output_root2 + '_events.fits'):
            raise Exception(
                'That events.fits file name is taken! Either delete the .fits '
                'file, or pick a new name.')
        if os.path.isfile(output_root2 + '_blends.fits'):
            raise Exception(
                'That blends.fits file name is taken! Either delete the .fits '
                'file, or pick a new name.')

    # Error handling/complaining if input types are not right.
    if type(output_root2) != str:
        raise Exception('output_root2 must be a string.')

    if type(radius_cut) != int:
        if type(radius_cut) != float:
            raise Exception('radius_cut must be an integer or a float.')

    if type(obs_time) != int:
        if type(obs_time) != float:
            raise Exception('obs_time must be an integer or a float.')

    if type(n_obs) != int:
        raise Exception('n_obs must be an integer.')

    if type(theta_frac) != int:
        if type(theta_frac) != float:
            raise Exception('theta_frac must be an integer or a float.')

    if seed is not None:
        if type(seed) != int:
            raise Exception('seed must be an integer.')

    ##########
    # Start of code
    #########

    # Set random seed
    np.random.seed(seed)

    t0 = time.time()

    # Initialize events_tmp and blends_tmp.
    events_tmp = None
    blends_tmp = None

    # Get the l and b from the HDF5 file.
    hf = h5py.File(hdf5_file, 'r')
    l_array = np.array(hf['long_bin_edges'])
    b_array = np.array(hf['lat_bin_edges'])
    hf.close()

    # Converts radius_cut from arcseconds into milliarcseconds
    radius_cut *= 1000.0

    # Set up the multiprocessing
    pool = Pool(n_proc)

    # Set up inputs to be able to be read by pool.map
    nll = len(l_array[:]) - 2
    nbb = len(b_array[:]) - 2

    llbb = itertools.product(range(nll), range(nbb))

    reps = nll * nbb

    hd = itertools.repeat(hdf5_file, reps)
    ot = itertools.repeat(obs_time, reps)
    no = itertools.repeat(n_obs, reps)
    rc = itertools.repeat(radius_cut, reps)
    tf = itertools.repeat(theta_frac, reps)
    br = itertools.repeat(blend_rad, reps)

    inputs = zip(llbb, hd, ot, no, rc, tf, br)

    ##########
    # Loop through galactic latitude and longitude bins. For each bin vertex,
    # take the nearest 4 bin samples and calculate microlensing events.
    # We do this to properly handle bin edges
    # (i.e. a sliding window analysis of 2x2 bins).
    # Duplicate events are removed.
    ##########
    # Should I use starmap_async?
    results = pool.starmap(_calc_event_time_loop, inputs)

    pool.close()
    pool.join()

    # Remove all the None values
    # (occurs for patches with less than 10 objects)
    results = [i for i in results if i is not None]

    results_ev = []
    results_bl = []

    for ii in range(len(results)):
        if results[ii] is not None:
            if results[ii][0] is not None:
                results_ev.append(results[ii][0])
            if results[ii][1] is not None:
                results_bl.append(results[ii][1])

    if len(results_ev) == 0:
        print('No events!')
        return
    else:
        events_tmp = np.concatenate(results_ev, axis=0)
        if len(results_bl) == 0:
            blends_tmp = np.array([])
        else:
            blends_tmp = np.concatenate(results_bl, axis=0)

    # Convert the events numpy recarray into an
    # Astropy Table for easier consumption.
    events_tmp = unique_events(events_tmp)
    events_final = Table(events_tmp)

    if len(results_bl) != 0:
        blends_tmp = unique_blends(blends_tmp)
    blends_final = Table(blends_tmp)

    # Save out file
    events_final.write(output_root2 + '_events.fits', overwrite=overwrite)
    blends_final.write(output_root2 + '_blends.fits', overwrite=overwrite)

    t1 = time.time()

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    radius_cut = radius_cut / 1000.0  # back to arcsec
    microlens_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    microlens_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=microlens_path).decode('ascii').strip()
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'
    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'hdf5_file , ' + hdf5_file + '\n'
    line2 = 'output_root2 , ' + output_root2 + '\n'
    line3 = 'radius_cut , ' + str(radius_cut) + ' , (arcsec)' + '\n'
    line4 = 'obs_time , ' + str(obs_time) + ' , (days)' + '\n'
    line5 = 'n_obs , ' + str(n_obs) + '\n'
    line6 = 'theta_frac , ' + str(theta_frac) + ' , (thetaE)' + '\n'
    line7 = 'blend_rad , ' + str(blend_rad) + ' , (arcsec)' + '\n'
    line8 = 'n_proc , ' + str(n_proc) + '\n'
    line9 = 'VERSION INFORMATION' + '\n'
    line10 = str(now) + ' : creation date' + '\n'
    line11 = microlens_hash + ' : microlens commit' + '\n'

    line12 = 'OTHER INFORMATION' + '\n'
    line13 = str(t1 - t0) + ' : total runtime (s)' + '\n'

    line14 = 'FILES CREATED' + '\n'
    line15 = output_root2 + '_events.fits : events file' + '\n'
    line16 = output_root2 + '_blends.fits : blends file' + '\n'

    with open(output_root2 + '_calc_events.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3,
                        line4, line5, line6, line7, line8, empty_line,
                        line9, dash_line, line10, line11, empty_line,
                        line12, dash_line, line13, empty_line, line14,
                        dash_line, line15, line16])

    print('Total runtime: {0:f} s'.format(t1 - t0))

    return


def _calc_event_time_loop(llbb, hdf5_file, obs_time, n_obs, radius_cut,
                          theta_frac, blend_rad):
    """
    Parameters
    ----------
    llbb : (int, int)
        Indices of (l,b) bin.

    obs_time, n_obs, radius_cut, theta_frac, blend_rad
    are all parameters of calc_events()

    Return
    ------
    events_llbb : array
        Array of the unique events for the particular (l,b) patch.

    blends_llbb : array
        Array of the unique blends for the particular (l,b) patch.

    """
    ####################
    # Loop through different time steps and figure out separations between
    # all possible pairs of stars. Trim down to "events", which consist of
    # those pairs that approach within one <radius_cut> of each other.
    # These will be the events we consider as candidate microlensing events.
    ####################

    # Initialize events_llbb and blends_llbb.
    events_llbb = None
    blends_llbb = None

    ll = llbb[0]
    bb = llbb[1]

    print('Working on loop ll, bb = ', ll, bb)
    name00 = 'l' + str(ll) + 'b' + str(bb)
    name01 = 'l' + str(ll) + 'b' + str(bb + 1)
    name10 = 'l' + str(ll + 1) + 'b' + str(bb)
    name11 = 'l' + str(ll + 1) + 'b' + str(bb + 1)

    hf = h5py.File(hdf5_file, 'r')
    bigpatch = np.hstack((hf[name00], hf[name01], hf[name10], hf[name11]))
    hf.close()

    # Skip patches with less than 10 objects
    if len(bigpatch) < 10:
        # continue
        return

    time_array = np.linspace(-1 * obs_time / 2.0, obs_time / 2.0, n_obs)

    for i in np.arange(len(time_array)):
        # Find potential lenses and sources that fall within radius cut.
        lens_id, sorc_id, r_t, sep, event_id1, c = _calc_event_cands_radius(bigpatch,
                                                                            time_array[i],
                                                                            radius_cut)

        # Calculate einstein radius and lens-source separation
        theta_E = einstein_radius(bigpatch['mass'][lens_id],
                                  r_t[lens_id], r_t[sorc_id])  # mas
        u = sep[event_id1] / theta_E

        # Trim down to those microlensing events that really get close enough
        # to hope that we can detect them. Trim on a Theta_E criteria.
        event_lbt = _calc_event_cands_thetaE(bigpatch, theta_E, u, theta_frac,
                                             lens_id, sorc_id, time_array[i])

        if event_lbt is not None:
            # Concatenate the current event table
            # (at this l, b, time) with the rest.
            if events_llbb is not None:
                events_llbb = np.hstack((events_llbb, event_lbt))
            else:
                events_llbb = event_lbt

            # Keep only unique events within our different time stamps
            events_llbb = unique_events(events_llbb)

            #########
            # Get blending.
            # Note 1: We are centering on the lens.
            # Note 2: We don't want to include the lens itself,
            # or the source, in the table.
            ##########
            blends_lbt = _calc_blends(bigpatch, c, event_lbt, blend_rad)

            if blends_lbt is not None:
                # Concatenate the current blend table (at this l, b, time)
                # with the rest.
                if blends_llbb is not None:
                    blends_llbb = np.hstack((blends_llbb, blends_lbt))
                else:
                    blends_llbb = blends_lbt

                # Keep only unique blends within our different time stamps
                blends_llbb = unique_blends(blends_llbb)

        # END of time loop

    return events_llbb, blends_llbb


def _calc_event_cands_radius(bigpatch, timei, radius_cut):
    """
    Get sources and lenses that pass the radius cut.

    Parameters
    ----------
    bigpatch : array
        Compilation of 4 .h5 datasets containing stars.

    timei : float
        Time at which to evaluate.

    radius_cut : float
        Parameter of calc_events().
        Converted to mas

    Return
    ------
    lens_id : array
        Indices into bigpatch that indicate lenses

    sorc_id : array
        Indices into bigpatch that indicate sources

    r_t : array
        Radial coordinates for all objects in bigpatch at time t

    sep : array
        Separation between lens-source pairs (mas)

    event_id1 : array
        Lens-source pairs where sep < radius_cut

    c : SkyCoord object
        Coordinates of all the stars.
    """
    # Propagate r, b, l positions forward in time.
    r_t = bigpatch['rad'] + timei * bigpatch['vr'] * kms_to_kpcday  # kpc
    b_t = bigpatch['glat'] + timei * bigpatch['mu_b'] * masyr_to_degday  # deg
    l_t = bigpatch['glon'] + timei * (bigpatch['mu_lcosb'] / np.cos(np.radians(bigpatch['glat']))) * masyr_to_degday  # deg

    ##########
    # Determine nearest neighbor in spherical coordinates.
    ##########
    c = SkyCoord(frame='galactic', l=l_t * units.deg, b=b_t * units.deg)

    # NOTE: dist has no actual meaning since
    # we didn't input distances from Earth.
    # It's an auto output. just ignore it.
    idx, sep, dist = coord.match_coordinates_sky(c, c, nthneighbor=2)

    # Converts separations to milliarcseconds
    sep = (sep.to(units.mas)) / units.mas

    ##########
    # Error checking: calculate how many duplicate (l, b) pairs there are.
    # (This is a problem for nearest neighbors.)
    ##########
    uni = np.unique((l_t, b_t), axis=1).shape[1]
    tot = len(l_t)
    dup = tot - uni
    if dup != 0:
        print('****************** WARNING!!! ********************')
        print('There are ' + str(dup) + ' duplicate (l, b) pairs.')
        print('**************************************************')

    ##########
    # Find all objects with a nearest neighbor within an angular distance
    # equal to sep. The index of the object and its nearest neighbor are
    # event_id1 and event_id2. (Indices correspond to those of idx and sep.)
    ##########
    # NOTE: event_id1/2 are indices into bigpatch
    event_id1 = np.where(sep < radius_cut)[0]
    event_id2 = idx[event_id1]

    ##########
    # We've got neighbors... figure out who's the lens and who's the source.
    ##########
    # NOTE: lens_id and sorc_id are indices into bigpatch
    idx_l1 = np.where(r_t[event_id1] < r_t[event_id2])[0]
    idx_l2 = np.where(r_t[event_id1] > r_t[event_id2])[0]

    lens_id = np.zeros(len(event_id1), dtype='int')
    sorc_id = np.zeros(len(event_id1), dtype='int')

    lens_id[idx_l1] = event_id1[idx_l1]
    sorc_id[idx_l1] = event_id2[idx_l1]

    lens_id[idx_l2] = event_id2[idx_l2]
    sorc_id[idx_l2] = event_id1[idx_l2]

    return lens_id, sorc_id, r_t, sep, event_id1, c


def _calc_event_cands_thetaE(bigpatch, theta_E, u, theta_frac, lens_id,
                             sorc_id, timei):
    """
    Get sources and lenses that pass the radius cut.

    Parameters
    ----------
    bigpatch : array
        Compilation of 4 .h5 datasets containing stars.

    theta_E : array
        Einstein radius of the events that pass the radius cut

    u : array
        Impact parameters at time t of the events that pass the radius cut

    theta_frac : float
        Parameter of calc_events()
        Another cut, in multiples of Einstein radii.

    lens_id : array
        Indices into bigpatch that indicate lenses

    sorc_id : array
        Indices into bigpatch that indicate sources

    timei : float
        Time at which to evaluate.

    Return
    ------
    event_lbt : array
        Lenses and sources at a particular time t.

    """
    # NOTE: adx is an index into lens_id or event_id (NOT bigpatch)
    adx = np.where(u < theta_frac)[0]
    if len(adx > 0):
        # Narrow down to unique pairs of stars... don't double calculate
        # an event.
        lens_sorc_id_pairs = np.stack((bigpatch['obj_id'][lens_id][adx],
                                       bigpatch['obj_id'][sorc_id][adx]),
                                      axis=-1)

        unique_returns = np.unique(lens_sorc_id_pairs,
                                   return_index=True, return_inverse=True,
                                   return_counts=True, axis=0)
        unique_tab = unique_returns[0]
        unique_indices = unique_returns[1]
        unique_inverse = unique_returns[2]
        unique_counts = unique_returns[3]

        # Define all tables we will need to calculate event properties.
        # These will all have a length of N_events.
        lens_table = bigpatch[lens_id][adx][unique_indices]
        sorc_table = bigpatch[sorc_id][adx][unique_indices]
        theta_E = theta_E[adx][unique_indices]
        u = u[adx][unique_indices]

        mu_b_rel = sorc_table['mu_b'] - lens_table['mu_b']  # mas/yr
        mu_lcosb_rel = sorc_table['mu_lcosb'] - lens_table['mu_lcosb']  # mas/yr
        mu_rel = np.sqrt(mu_b_rel ** 2 + mu_lcosb_rel ** 2)  # mas/yr
        t_event = np.ones(len(mu_rel), dtype=float) * timei  # days

        # This is all the events for this l, b, time
        # Loop through the lens table and append '_L' to the end of each field
        lens_rename_dct = {}
        for name in lens_table.dtype.names:
            lens_rename_dct[name] = name + '_L'
        lens_table = rfn.rename_fields(lens_table, lens_rename_dct)

        # Loop through the source table and append '_S' to the end of each field
        sorc_rename_dct = {}
        for name in sorc_table.dtype.names:
            sorc_rename_dct[name] = name + '_S'
        sorc_table = rfn.rename_fields(sorc_table, sorc_rename_dct)

        # Combine the lens and source tables into the events table
        event_lbt = rfn.merge_arrays((lens_table, sorc_table),
                                     flatten=True)

        # Add additional microlensing parameters to the events table
        event_lbt = rfn.append_fields(event_lbt, 'theta_E',
                                      theta_E, usemask=False)
        event_lbt = rfn.append_fields(event_lbt, 'u0',
                                      u, usemask=False)
        event_lbt = rfn.append_fields(event_lbt, 'mu_rel',
                                      mu_rel, usemask=False)
        event_lbt = rfn.append_fields(event_lbt, 't0',
                                      t_event, usemask=False)
        return event_lbt

    else:
        return None


def _calc_blends(bigpatch, c, event_lbt, blend_rad):
    """
    Create a table containing the blended stars for each event.
    Note 1: We are centering on the lens.
    Note 2: We don't want to include the lens itself,
    or the source, in the table.

    Parameters
    ----------
    bigpatch : array
        Compilation of 4 .h5 datasets containing stars.

    c : SkyCoord object
        Coordinates of all the stars.

    event_lbt : array
        Lenses and sources at a particular time t.

    blend_rad : float
        Parameter of calc_events().

    Return
    ------
    blends_lbt : array
        Array of neighbor stars for each lens-source pair.

    """
    ##########
    # Get the cached KD-Tree to make things run faster.
    ##########
    # This way, we don't have to remake a tree that already exists.
    # (We didn't just do the neighbor calculation initially, because
    # it would be expensive to hold onto all the unnecessary neighbors)
    kdtree_cache = c.cache['kdtree_sky']

    # Define the center of the blending disk (the lens)
    coords_lbt = SkyCoord(frame='galactic',
                          l=np.array(event_lbt['glon_L']) * units.deg,
                          b=np.array(event_lbt['glat_L']) * units.deg)

    ##########
    # Replicate astropy's search_around_sky.
    ##########
    # Make the coordinates to query around
    seplimit = blend_rad * units.arcsec
    coords1 = coords_lbt
    coords1 = coords1.transform_to(c)
    urepr1 = coords1.data.represent_as(UnitSphericalRepresentation)
    ucoords1 = coords1.realize_frame(urepr1)
    cartxyz1 = ucoords1.cartesian.xyz
    flatxyz1 = cartxyz1.reshape((3, np.prod(cartxyz1.shape) // 3))

    # Define the query distance.
    r_kdt = (2 * np.sin(Angle(seplimit) / 2.0)).value

    # Query ball against the existing (cached) tree.
    # NOTE: results is an array of lists.
    results = kdtree_cache.query_ball_point(flatxyz1.T, r_kdt)

    # Figure out the number of blends for each lens.
    blend_lens_obj_id = []
    blend_sorc_obj_id = []
    blend_neigh_obj_id = []
    blend_neigh_idx = []
    sep_LN_list = []

    for ii in range(len(results)):
        # results indexes into bigpatch.
        # ii corresponds to coords_lbt.
        if len(results[ii]) == 0:
            continue

        # bidx indexes into results.
        # It should be that len(results[ii]) == len(bidx) + 2 (we get rid of source and lens.)
        # neighbor star object id
        nid = bigpatch['obj_id'][results[ii]]
        # lens star object id
        lid = np.array(event_lbt['obj_id_L'][ii])
        # source star object id
        sid = np.array(event_lbt['obj_id_S'][ii])

        # Fetch the things that are not the lens and the source.
        bidx = np.where((nid != lid) &
                        (nid != sid))[0]

        if len(bidx) == 0:
            continue

        # Make a list of the lens and source IDs for each neighbor... these are just repeats.
        tmp_lens_id = np.repeat(lid, len(bidx)).tolist()
        tmp_sorc_id = np.repeat(sid, len(bidx)).tolist()

        # Calculate the distance from lens to each neighbor.
        lens_lb = SkyCoord(frame='galactic',
                           l=np.array(event_lbt['glon_L'][ii]) * units.deg,
                           b=np.array(event_lbt['glat_L'][ii]) * units.deg)
        neigh_lb = SkyCoord(frame='galactic',
                            l=bigpatch['glon'][results[ii]][bidx] * units.deg,
                            b=bigpatch['glat'][results[ii]][bidx] * units.deg)
        sep_LN = lens_lb.separation(neigh_lb)
        sep_LN = (sep_LN.to(units.arcsec)) / units.arcsec

        # Add the non-lens, non-source blended object IDs to a list.
        # Add the lens and the source object ID to a new list as well.
        # Add the index of the neighbors.
        blend_neigh_obj_id.extend(nid[bidx].tolist())
        blend_lens_obj_id.extend(tmp_lens_id)
        blend_sorc_obj_id.extend(tmp_sorc_id)
        blend_neigh_idx.extend([results[ii][bb] for bb in bidx])
        sep_LN_list.extend(sep_LN.value.tolist())

    # Convert our lists into arrays.
    blend_neigh_obj_id = np.array(blend_neigh_obj_id)
    blend_lens_obj_id = np.array(blend_lens_obj_id)
    blend_sorc_obj_id = np.array(blend_sorc_obj_id)
    blend_neigh_idx = np.array(blend_neigh_idx)
    sep_LN_list = np.array(sep_LN_list)

    if len(blend_neigh_obj_id) > 0:
        # Grab the rows of bigpatch that are neighbors
        blends_lbt = bigpatch[blend_neigh_idx]

        # Append '_N' to each column in blends_lbt
        blends_rename_dct = {}
        for name in blends_lbt.dtype.names:
            blends_rename_dct[name] = name + '_N'
        blends_lbt = rfn.rename_fields(blends_lbt, blends_rename_dct)

        # Add additional columns into blends_lbt
        blends_lbt = rfn.append_fields(blends_lbt, 'obj_id_L',
                                       blend_lens_obj_id, usemask=False)
        blends_lbt = rfn.append_fields(blends_lbt, 'obj_id_S',
                                       blend_sorc_obj_id, usemask=False)
        blends_lbt = rfn.append_fields(blends_lbt, 'sep_LN',
                                       sep_LN_list, usemask=False)
    else:
        blends_lbt = None

    return blends_lbt


def unique_events(event_table):
    """
    Given an event table, there might be a single microlensing event listed
    multiple times (it might have been observed at different timesteps, or
    the nearest neighbor pair might have been symmetric for source and lens.)
    This function will eliminate duplicates, only keeping an event once. It
    is picked to be the particular observed event with the smallest source-
    lens separation.

    Parameters
    ---------
    event_table : numpy array 
        A table with all the events. There are equal numbers of columns
        containing info about the source and the lens, and four additional
        columns with theta_E, u0, mu_rel, and t0. The number of rows
        corresponds to the number of events.

    Return
    ------
    new_event_table : numpy array
        Same as event_table, but all duplicate events have been trimmed out,
        such that each event only is listed once (at the observed time where
        the source-lens separation is smallest.)

    """
    # Pull the unique ID numbers for the lens and source and put them into a
    # table of 2 x N_lenses.
    lens_uid = event_table['obj_id_L']
    sorc_uid = event_table['obj_id_S']
    events_uid = np.swapaxes(np.vstack((lens_uid, sorc_uid)), 0, 1)

    # Determine if we have unique events (and how many duplicates there are).
    unique_returns = np.unique(events_uid,
                               return_index=True, return_inverse=True,
                               return_counts=True, axis=0)
    unique_tab = unique_returns[0]
    unique_indices = unique_returns[1]
    unique_inverse = unique_returns[2]
    unique_counts = unique_returns[3]

    new_event_table = event_table[unique_indices]

    # Check for duplicate events and keep the one with the closest u.
    dpdx = np.where(unique_counts > 1)[0]
    for ii in range(len(dpdx)):
        # Fetch the duplicates for this event.
        dup_idx = np.where(unique_inverse == dpdx[ii])[0]
        dup_events = event_table[dup_idx]
        min_idx = np.argmin(dup_events['u0'])
        new_event_table[dpdx[ii]] = event_table[dup_idx[min_idx]]

    return new_event_table


def unique_blends(blend_table):
    """
    Given an blends table, there might be a single lens-source-neighbor triple
    multiple times (it might have been observed at different timesteps.)
    This function will eliminate duplicates, only keeping an event once. It
    is picked to be the first occurence.

    Parameters
    ---------
    blend_table : blend array 
        A table with all the events. There is 1 column with the unique
        source ID, 1 with the unique lens ID lens, 1 with the lens-neighbor
        separation, and the remaining columns contain info about the neighbors.


    Return
    ------
    new_blend_table : numpy array
        Same as blend_table, but all duplicate events have been trimmed out,
        such that each event only is listed once (at the observed time where
        the source-lens separation is smallest.)

    """
    # Pull the unique ID numbers for the lens, source, and neighbors and put
    # them into a table of 3 x N_blends.
    lens_obj_id = blend_table['obj_id_L']
    sorc_obj_id = blend_table['obj_id_S']
    neigh_obj_id = blend_table['obj_id_N']

    triples_obj_id = np.swapaxes(np.vstack((lens_obj_id,
                                            sorc_obj_id,
                                            neigh_obj_id)),
                                 0, 1)

    # Determine if we have unique events (and how many duplicates there are).
    # We will keep the first occurence.
    unique_returns = np.unique(triples_obj_id, return_index=True, axis=0)
    unique_tab = unique_returns[0]
    unique_indices = unique_returns[1]

    unique_blend_tab = blend_table[unique_indices]

    return unique_blend_tab


def calc_diff_limit_blend(event_fits_file, blend_fits_file, blend_rad):
    """
    Given a table of events and blends, calculate what the blending would be
    for a smaller blending radius (e.g. table is for seeing limited,
    and you want what the diffraction limit blend is.)
    """
    diff_limit_blend = None

    event = Table.read(event_fits_file)
    blend = Table.read(blend_fits_file)

    for i in range(len(blend)):
        lens_id = blend[i]['obj_id_L']
        sorc_id = blend[i]['obj_id_S']
        event_idx = np.where(
            (event['obj_id_L'] == lens_id) & (event['obj_id_S'] == sorc_id))[0]

        l_L = event[event_idx]['glon_L']
        b_L = event[event_idx]['glat_L']
        c_L = SkyCoord(frame='galactic', l=l_L * units.deg, b=b_L * units.deg)

        l_N = blend[i]['glon_N']
        b_N = blend[i]['glat_N']
        c_N = SkyCoord(frame='galactic', l=l_N * units.deg, b=b_N * units.deg)

        sep = (c_L.separation(c_N)).arcsec

        if sep < blend_rad:
            if diff_limit_blend is not None:
                diff_limit_blend = vstack((diff_limit_blend, blend[i]))
            else:
                diff_limit_blend = blend[i]

    diff_limit_blend.write('diff_limit_blend_' + str(blend_rad) + '.fits')

    return


def reduce_blend_rad(blend_tab, new_blend_rad, output_root, overwrite=False):
    """
    Creates a new blend table for some blending radius r_new
    that is smaller than the original blend radius r_orig,
    i.e. r_new < r_orig. Also makes a corresponding log and events.

    Parameters
    ----------
    blend_tab : str
        The name of the blend table.

    new_blend_rad : float or int
        The new (smaller) blend radius.
        Units are in ARCSECONDS.

    output_root : str
        The name for the new blend table
        (and corresponding event table)

    Return
    ------
    new_blend : .fits table
        New table with smaller blend radius.

    """
    input_root = blend_tab.rstrip('_blends.fits')

    event_tab_name = input_root + '_events.fits'
    event_log_name = input_root + '_calc_events.log'

    new_event_tab_name = output_root + '_events.fits'
    new_event_log_name = output_root + '_calc_events.log'
    new_blend_tab_name = output_root + '_blends.fits'

    os.system('cp %s %s' % (event_tab_name, new_event_tab_name))
    os.system('cp %s %s' % (event_log_name, new_event_log_name))

    now = datetime.datetime.now()
    dash_line = '-----------------------------' + '\n'

    line0 = 'NEW BLEND RADIUS INFO AND FILES CREATED' + '\n'
    line1 = 'new blend rad : ' + str(new_blend_rad) + ' (arcsec)' + '\n'
    line2 = 'blend file : ' + new_blend_tab_name + '\n'
    line3 = 'event file : ' + new_event_tab_name + '\n'
    line4 = 'creation time: ' + str(now) + '\n'

    with open(new_event_log_name, 'a') as my_file:
        my_file.writelines(['\n', line0, dash_line,
                            line1, line2, line3, line4])

    old_blends = Table.read(blend_tab)
    good_idx = np.where(old_blends['sep_LN'] < new_blend_rad)[0]
    new_blends = old_blends[good_idx]

    new_blends.write(new_blend_tab_name, overwrite=overwrite)

    return


############################################################################
######### Refined event rate calculation and associated functions ##########
############################################################################


def _convert_photometric_99_to_nan(table, photometric_system='ubv'):
    for name in table.colnames:
        if ('exbv' in name) or (photometric_system in name):
            cond = np.where(table[name] == -99)[0]
            table[name][cond] = np.nan


def refine_events(input_root, filter_name, photometric_system, red_law,
                  overwrite=False,
                  output_file='default'):
    """
    Takes the output Astropy table from calc_events, and from that
    calculates the time of closest approach. Will also return source-lens
    separation at this time.

    Parameters
    ----------
    input_root : str
        The root path and name of the *_events.fits and *_blends.fits.
        Don't include those suffixes yet.

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    photometric_system : str
        The name of the photometric system in which the filter exists.

    red_law : str
        The name of the reddening law to use from PopStar.

    Optional Parameters
    -------------------
    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    Output:
    ----------
    A file will be created named
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits that contains all the
    same objects, only now with lots of extra
    columns of data.

    """
    ##########
    # Error handling: check whether files exist and
    # whether input types are correct.
    ##########

    # Error handling/complaining if input types are not right.
    if type(input_root) != str:
        raise Exception('input_root must be a string.')

    if type(filter_name) != str:
        raise Exception('filter_name must be a string.')

    if type(red_law) != str:
        raise Exception('red_law must be a string.')

    if type(output_file) != str:
        raise Exception('output_file must be a string.')

    # Check to see that the filter name, photometric system, red_law are valid
    if photometric_system not in photometric_system_dict:
        exception_str = 'photometric_system must be a key in ' \
                        'photometric_system_dict. ' \
                        'Acceptable values are : '
        for photometric_system in photometric_system_dict:
            exception_str += '%s, ' % photometric_system
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    if filter_name not in photometric_system_dict[photometric_system]:
        exception_str = 'filter_name must be a value in ' \
                        'photometric_system_dict[%s]. ' \
                        'Acceptable values are : ' % photometric_system
        for filter_name in photometric_system_dict[photometric_system]:
            exception_str += '%s, ' % filter_name
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    key = photometric_system + '_' + filter_name
    if red_law not in filt_dict[key]:
        raise Exception('%s not in filt_dict[%s]' % (red_law, key))

    # Check if .fits file exists already. If it does, throw an error message
    # to complain and exit.
    # Error handling.
    if not overwrite and os.path.isfile(output_file):
        raise Exception('That refined_events.fits file name is taken! '
                        'Either delete the .fits file, or pick a new name.')

    t_0 = time.time()

    if output_file == 'default':
        output_file = '{0:s}_refined_events_{1:s}_{2:s}_{3:s}.fits'.format(input_root,
                                                                           photometric_system,
                                                                           filter_name,
                                                                           red_law)

    event_fits_file = input_root + '_events.fits'
    blend_fits_file = input_root + '_blends.fits'
    log_file = input_root + '_calc_events.log'

    event_tab = Table.read(event_fits_file)
    blend_tab = Table.read(blend_fits_file)

    # If photometric fields contain -99, convert to nan
    _convert_photometric_99_to_nan(event_tab, photometric_system)
    _convert_photometric_99_to_nan(blend_tab, photometric_system)

    # Only keep events with luminous sources
    event_tab = event_tab[~np.isnan(event_tab[photometric_system + '_' + filter_name + '_S'])]

    with open(log_file) as my_file:
        for num, line in enumerate(my_file):
            if 'obs_time' in line:
                obs_time = line.split(',')[1]
                obs_time = int(obs_time)
                break

    # Calculate time and separation at closest approach, add to table
    # NOTE: calc_closest_approach modifies the event table!
    # It trims out events that peak outside the survey range!
    print('Original candidate events: ', len(event_tab))
    u0, t0 = calc_closest_approach(event_tab, obs_time)
    print('Candidate events in survey window: ', len(event_tab))
    event_tab['t0'] = t0  # days
    event_tab['u0'] = u0
    if len(event_tab) == 0:
        print('No events!')
        return

    ##########
    # Calculate apparent magnitudes
    ##########

    # Einstein crossing time
    # THIS HAS TO GO BEFORE _CALC_OBSERVABLES
    t_E = event_tab['theta_E'] / event_tab['mu_rel']  # yr
    t_E *= 365.25  # days
    event_tab['t_E'] = t_E  # days

    # Add stuff to event_tab... shouldn't have any direct outputs
    _calc_observables(filter_name, red_law, event_tab, blend_tab, photometric_system)

    # Relative parallax
    pi_rel = event_tab['rad_L'] ** -1 - event_tab['rad_S'] ** -1
    event_tab['pi_rel'] = pi_rel  # mas

    # Microlensing parallax
    pi_E = pi_rel / event_tab['theta_E']
    event_tab['pi_E'] = pi_E  # dim'less

    event_tab.write(output_file, overwrite=overwrite)

    t_1 = time.time()

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    microlens_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    microlens_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=microlens_path).decode('ascii').strip()
    popstar_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                           cwd=popstar_path).decode('ascii').strip()
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'input_root : ' + input_root + '\n'
    line2 = 'filter_name : ' + filter_name + '\n'
    line3 = 'red_law : ' + red_law + '\n'

    line4 = 'VERSION INFORMATION' + '\n'
    line5 = str(now) + ' : creation date' + '\n'
    line6 = popstar_hash + ' : PopStar commit' + '\n'
    line7 = microlens_hash + ' : microlens commit' + '\n'

    line8 = 'OTHER INFORMATION' + '\n'
    line9 = str(t_1 - t_0) + ' : total runtime (s)' + '\n'

    line10 = 'FILES CREATED' + '\n'
    line11 = output_file + ' : refined events'

    with open(input_root + '_refined_events_' + photometric_system + '_' + filter_name + '_' + red_law + '.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, empty_line,
                        line4, dash_line, line5, line6, line7, empty_line,
                        line8, dash_line, line9, empty_line,
                        line10, dash_line, line11])

    print('Total runtime: {0:f} s'.format(t_1 - t_0))
    return


def calc_closest_approach(event_tab, survey_duration):
    """
    Calculate the point of closest approach, during a given interval.
    Return the position and time of closest approach.
    This algorithm comes from a little writeup called
    "Point of Closest Approach for the Linear Motion of Two Particles in 2-D"

    Parameters
    ----------
    event_tab : astropy table
        This should take the astropy table that comes from the
        output of calc_events.

    survey_duration : int or float
        Survey duration, in days

    Return
    ------
    u0 : array (same length as event_tab)
        Minimum separation (normalized to Einstein radius)

    t0 : array (same length as event_tab)
        Time at minimum separation (days)

    """
    # Survey start and end time
    start_time = -1 * survey_duration / 2.0  # days
    end_time = survey_duration / 2.0  # days

    # Parametrization coefficients
    cos_glat_S = np.cos(np.radians(event_tab['glat_S']))
    cos_glat_L = np.cos(np.radians(event_tab['glat_L']))

    l_0 = (event_tab['glon_S'] - event_tab['glon_L']) * cos_glat_L * 3600.0 * 1.0e3  # position, mas
    l_1 = (event_tab['mu_lcosb_S'] - event_tab['mu_lcosb_L'])  # pm, mas/yr
    b_0 = (event_tab['glat_S'] - event_tab['glat_L']) * 3600.0 * 1.0e3  # position, mas
    b_1 = (event_tab['mu_b_S'] - event_tab['mu_b_L'])  # pm, mas/yr

    # Calculate the critical time (time of closest approach), global minimum
    t_crit = -(l_0 * l_1 + b_0 * b_1) / (l_1 ** 2 + b_1 ** 2)  # years
    t_crit *= 365.25  # days

    # Only include events that peak in survey window,
    # since we can't easily fit parameters for events where we don’t have the peak.
    # FIXME: more efficient way to do this??
    good = np.where((t_crit > start_time) & (t_crit < end_time))[0]
    bad = np.where((t_crit <= start_time) | (t_crit >= end_time))[0]

    event_tab.remove_rows(bad)

    # Initialize the time of closest approach vector, minimum of survey
    # NOTE: THIS WILL OVERWRITE CURRENT VALUES
    t0 = np.zeros(len(event_tab))
    u0 = np.zeros(len(event_tab))

    # Assign values to closest approach values
    t0 = t_crit[good]
    u0 = calc_distance(event_tab, t0)

    return u0, t0


def calc_distance(event_tab, time):
    """
    Calculate the separation of two different objects at some given time.

    Parameters
    ----------
    event_tab : astropy table

    time : int/float

    Return
    ------
    u : array (same length as astropy table)

    """
    # Calculate sep
    cos_glat_S = np.cos(np.radians(event_tab['glat_S']))
    cos_glat_L = np.cos(np.radians(event_tab['glat_L']))

    l_L = (event_tab['glon_L'] + time * event_tab['mu_lcosb_L'] * masyr_to_degday / cos_glat_L)
    b_L = (event_tab['glat_L'] + time * event_tab['mu_b_L'] * masyr_to_degday)
    l_S = (event_tab['glon_S'] + time * event_tab['mu_lcosb_S'] * masyr_to_degday / cos_glat_S)
    b_S = (event_tab['glat_S'] + time * event_tab['mu_b_S'] * masyr_to_degday)

    c_L = SkyCoord(frame='galactic', l=l_L * units.deg, b=b_L * units.deg)
    c_S = SkyCoord(frame='galactic', l=l_S * units.deg, b=b_S * units.deg)

    sep = (c_L.separation(c_S)).mas
    u = sep / event_tab['theta_E']

    return u


def calc_blend_and_centroid(filter_name, red_law, blend_tab, photometric_system='ubv'):
    """
    Given the absolute magnitudes of a bunch of things,
    calculate their blended apparent magnitude and flux.
    Also calculate the centroid of these things.

    Filter name is 'j', 'i', etc.
    red_law is Damineli16, Schlegel99, etc.
    """
    f_i = filt_dict[photometric_system + '_' + filter_name][red_law]

    # Calculate apparent magnitudes
    app_N = calc_app_mag(blend_tab['rad_N'],
                         blend_tab[photometric_system + '_' + filter_name + '_N'],
                         blend_tab['exbv_N'], f_i)

    # Convert absolute magnitudes to fluxes, and fix bad values
    flux_N = 10 ** (app_N / -2.5)
    flux_N = np.nan_to_num(flux_N)

    # Get total flux
    flux_N_tot = np.sum(flux_N)

    # Get centroid
    # If there is no flux, then the centroid is set to (l, b) = (0, 0).
    if flux_N_tot == 0:
        cent_l = 0.0
        cent_b = 0.0
    else:
        cent_l = np.sum(flux_N * blend_tab['glon_N']) / flux_N_tot
        cent_b = np.sum(flux_N * blend_tab['glat_N']) / flux_N_tot

    # Total blended magnitude
    app_blended = -2.5 * np.log10(flux_N_tot)

    return app_blended, flux_N_tot, cent_l, cent_b


def _calc_observables(filter_name, red_law, event_tab, blend_tab, photometric_system='ubv'):
    """
    Calculate a bunch of observable quantities we get out from microlensing

    Parameters
    ----------
    filter_name : str

    event_tab :

    blend_tab :

    Return
    ------
    """
    f_i = filt_dict[photometric_system + '_' + filter_name][red_law]

    # Calculate apparent magnitude of lens and source, and fix bad values
    app_S = calc_app_mag(event_tab['rad_S'],
                         event_tab[photometric_system + '_' + filter_name + '_S'],
                         event_tab['exbv_S'], f_i)
    app_L = calc_app_mag(event_tab['rad_L'],
                         event_tab[photometric_system + '_' + filter_name + '_L'],
                         event_tab['exbv_L'], f_i)
    event_tab[photometric_system + '_' + filter_name + '_app_S'] = app_S
    event_tab[photometric_system + '_' + filter_name + '_app_L'] = app_L

    # Convert absolute magnitude to fluxes, and fix bad values
    flux_L = 10 ** (app_L / -2.5)
    flux_S = 10 ** (app_S / -2.5)

    flux_L = np.nan_to_num(flux_L)
    flux_S = np.nan_to_num(flux_S)

    # Find the blends.
    LS_pairs = np.stack((event_tab['obj_id_L'], event_tab['obj_id_S']),
                        axis=-1)
    blend_pairs = np.stack((blend_tab['obj_id_L'], blend_tab['obj_id_S']),
                           axis=-1)

    flux_N = np.zeros(len(app_L))
    event_tab['cent_glon_' + filter_name + '_N'] = np.zeros(len(app_L))
    event_tab['cent_glat_' + filter_name + '_N'] = np.zeros(len(app_L))

    # Get the unique blend_pairs and the first instance of each.
    if len(blend_pairs) > 0:
        uni_blends, uni_bidx = np.unique(blend_pairs, return_index=True,
                                         axis=0)
        # Add a "dummy" idx onto uni_bidx so that I can do idx+1 later.
        uni_bidxx = np.append(uni_bidx, len(blend_pairs))
    else:
        # obj_id always >= 0, so these negative values serve as dummies that
        # cannot be matched in the below `where` statement
        uni_blends = np.array([[-1, -1]])
        uni_bidxx = np.array([])

    # FIXME: Put in a check that uni_bidx is monotonically increasing???
    # np.sort(uni_bidx) - uni_bidx == np.zeros(len(uni_bidx)) ???

    # FIXME: This seems more complicated than necessary...
    for pp in range(len(LS_pairs)):
        # Search for where the blend pairs have the nth unique LS value
        idx = np.where(uni_blends == LS_pairs[pp])[0]
        if len(idx) > 0:
            start = uni_bidxx[idx][0]
            end = uni_bidxx[idx + 1][0]
            app_blended, flux_N_tot, cent_l, cent_b = calc_blend_and_centroid(filter_name,
                                                                              red_law,
                                                                              blend_tab[start:end],
                                                                              photometric_system)
            flux_N[pp] = flux_N_tot
            event_tab['cent_glon_' + filter_name + '_N'][pp] = cent_l
            event_tab['cent_glat_' + filter_name + '_N'][pp] = cent_b

    # Total blended flux in i-band
    flux_tot = flux_L + flux_S + flux_N

    # Total blended magnitude, neighbor
    app_N = np.zeros(len(flux_N))
    bad_N_idx = np.where(flux_N == 0)[0]
    good_N_idx = np.where(flux_N != 0)[0]
    app_N[good_N_idx] = -2.5 * np.log10(flux_N[good_N_idx])
    app_N[bad_N_idx] = np.full(len(bad_N_idx), np.nan)

    event_tab[photometric_system + '_' + filter_name + '_app_N'] = app_N

    # Total blended magnitude
    app_LSN = -2.5 * np.log10(flux_tot)
    event_tab[photometric_system + '_' + filter_name + '_app_LSN'] = app_LSN

    # Bump amplitude (in magnitudes)
    delta_m = calc_bump_amp(event_tab['u0'], flux_S, flux_L, flux_N)
    event_tab['delta_m_' + filter_name] = delta_m

    # Calculate the blend fraction
    # (commonly used in microlensing): f_source / f_total
    f_blend = flux_S / flux_tot
    event_tab['f_blend_' + filter_name] = f_blend

    return


##################################################################
############ Reading/writing and format functions  ###############
##################################################################

def make_ebf_log(ebf_table):
    """
    Converts log from Galaxia ebf output into dictionary

    Parameters
    ----------
    ebf_table : processed ebf file
        The ebf file that galaxia outputs, AFTER it's been read

    Returns
    -------
    ebf_log : dictionary
        The ebf file log, in dictionary form

    """
    if type(ebf_table) == dict:
        log_list = ebf_table['log'][0].split(b'\n')

    if ebf_table[-4:] == '.ebf':
        log_list = ebf.read(ebf_table, '/log')[0].split(b'\n')

    if type(ebf_table) == np.bytes_:
        log_list = ebf_table.split(b'\n')

    ebf_log = {}
    for ii in range(len(log_list)):
        if len(log_list[ii]) == 0:
            continue

        if log_list[ii].startswith(b'#'):
            continue

        log_list_entry = log_list[ii].split()

        ebf_log[log_list_entry[0].decode('utf-8')] = log_list_entry[1].decode(
            'utf-8')

    return ebf_log


def make_label_file(h5file_name, overwrite=False):
    """
    Given an output root for an .h5 file, creates a table of
    dataset name, l, b, and number of objects.
    Writes out the Astropy table as a .fits file.

    Parameters
    ----------
    h5file_name : string
        The path and name of the input h5file containing the
        population synthesis results. We will read this in and
        make a new output file entitled:
        <h5file_name>_label.fits (after removing the *.h5).

    Return
    ------
    label_file : Astropy table
        Table containing the dataset name, and corresponding l, b, and number of objects.

    """
    data_dict = {'file_name': [], 'long_start': [], 'long_end': [],
                 'lat_start': [], 'lat_end': [], 'objects': [],
                 'N_stars': [], 'N_WD': [], 'N_NS': [], 'N_BH': []}

    hf = h5py.File(h5file_name + '.h5', 'r')
    l_array = hf['long_bin_edges']
    b_array = hf['lat_bin_edges']

    for ll in range(len(l_array) - 1):
        for bb in range(len(b_array) - 1):
            dset_name = 'l' + str(ll) + 'b' + str(bb)
            dataset = hf[dset_name]

            if len(dataset) == 0:
                N_stars = 0
                N_WD = 0
                N_NS = 0
                N_BH = 0
            else:
                N_stars = np.sum(dataset['rem_id'] == 0)
                N_WD = np.sum(dataset['rem_id'] == 101)
                N_NS = np.sum(dataset['rem_id'] == 102)
                N_BH = np.sum(dataset['rem_id'] == 103)

            data_dict['file_name'].append(dset_name)
            data_dict['long_start'].append(l_array[ll])
            data_dict['long_end'].append(l_array[ll + 1])
            data_dict['lat_start'].append(b_array[bb])
            data_dict['lat_end'].append(b_array[bb + 1])
            data_dict['objects'].append(dataset.shape[0])
            data_dict['N_stars'].append(N_stars)
            data_dict['N_WD'].append(N_WD)
            data_dict['N_NS'].append(N_NS)
            data_dict['N_BH'].append(N_BH)

    hf.close()

    label_file = Table(data_dict, names=('file_name', 'long_start', 'long_end',
                                         'lat_start', 'lat_end', 'objects',
                                         'N_stars', 'N_WD', 'N_NS', 'N_BH'))
    label_file['long_start'].format = '08.3f'
    label_file['long_end'].format = '08.3f'
    label_file['lat_start'].format = '07.3f'
    label_file['lat_end'].format = '07.3f'

    now = datetime.datetime.now()
    label_file.meta['label'] = 'label.fits file creation time: ' + str(now)

    label_file.write(h5file_name + '_label.fits', overwrite=overwrite)

    return

###########################################################################
############ General formulas, conversions, and calculations ##############
###########################################################################

def heliocentric_to_galactic(x, y, z):
    """
    Converts from heliocentric coordinates to galactic coordinates.

    Parameters
    ----------
    x, y, z : float or array
        Heliocentric coordinates x, y, and z (in kpc)

    Returns
    -------
    r, b, l : float or array
        Galactic coordinates r, b, and l (in kpc and degrees)

    """
    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5

    # np.arccos will return a value between 0 and pi
    # so b will be between -90 and 90 degrees
    b = -np.degrees(np.arccos(z / r)) + 90

    # np.arctan2 will return a value between -pi and pi
    # arctan2 takes care of choosing the correct branch of the arctan function
    # %360 will return it between 0 and 360 degrees
    l = np.degrees((np.arctan2(y, x))) % 360

    return r, b, l


def galactic_to_heliocentric(r, b, l):
    """
    Converts from galactic coordinates to heliocentric coordinates.

    Parameters
    ---------
    r, b, l : float or array
        Galactic coordinates r, b and l (in kpc and degrees)

    Returns
    -------
    x, y, z : float or array
        Heliocentric coordinates x, y, and z (in kpc)

    """
    # Convert input from degrees to radians.
    l = np.radians(l)
    b = np.radians(b)

    x = r * np.cos(b) * np.cos(l)
    y = r * np.cos(b) * np.sin(l)
    z = r * np.sin(b)

    return x, y, z


def einstein_radius(M, d_L, d_S):
    """
    Calculates the einstein radius, in mas

    Parameters
    ----------
    M : float
        Lens mass, in solar masses

    d_L : float
        Distance to lens, in kpc

    d_S : float
        Distance to source, in kpc

    Return
    ------
        Einstein radius, in mas
    """
    return 2.85 * M ** 0.5 * (1 / d_L - 1 / d_S) ** 0.5


def calc_sph_motion(vx, vy, vz, r, b, l):
    """
    Calculate velocities in the r directions and proper motions
    in l, b directions.

    Parameters
    ----------
    vx, vy, vz : float or array
        Heliocentric velocities vx, vy, and vz (in km/s)

    r, b, l : float or array
        Galactic coordinates r, b, and l (in kpc and degrees)

    Return
    ------
    vr, mu_b, mu_lcosb : float or array (in km/s and mas/yr)
        Radial velocity and proper motions

    """
    b = np.radians(b)
    l = np.radians(l)

    cosb = np.cos(b)
    sinb = np.sin(b)
    cosl = np.cos(l)
    sinl = np.sin(l)

    vr = vz * sinb + cosb * (vx * cosl + vy * sinl)
    vb = vz * cosb - sinb * (vx * cosl + vy * sinl)
    vl = vy * cosl - vx * sinl

    mu_b = vb / r
    mu_lcosb = vl / r

    ##########
    # Convert between radians*(km/s)/kpc into mas/year.
    # 1 kpc = 3.086 * 10^16 km, 1 year = 3.154 * 10^7 s, 2.063 * 10^8 mas = 1 rad
    # 1 radian*(km/s)/kpc = 0.2108 mas/yr
    #########
    conversion_fact = 0.2108
    mu_b = mu_b * conversion_fact
    mu_lcosb = mu_lcosb * conversion_fact

    return vr, mu_b, mu_lcosb


def calc_normalized_counts(mag):
    """
    CHANGE THIS CODE IN THE FUTURE TO TAKE IN DIFFERENT ZERO POINTS!
    Right now this is only applicable for OGLE I band
    See confluence wiki for where these values come from...
    """
    a = -0.67815949
    b = 15.37993393
    counts = 10 ** (b + a * mag)
    return counts


def calc_magnification(u):
    """
    Calculate the magnification factor A(u)

    Parameters
    ----------
    u : float or array
        Dimensionless lens-source separation, normalized in Einstein radius units

    Return
    ------
    Magnification factor : float or array

    """
    return (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))


def calc_delta_c(u, thetaE):
    """
    Calculate the maximum centroid shift for a dark lens,
    no neighbors
    """
    delta_c = u * thetaE / (u ** 2 + 2)

    return delta_c


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

    Return
    ------
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

    Return
    ------
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

    Return
    ------
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

    Return
    ------
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


def calc_bump_amp(u0, f_S, f_L, f_N):
    """
    Calculate the "bump" amplitude, given the minimum separation and the fluxes
    of the (unmagnified) source, lens, and neighbors.
    The bump amplitude (mags) is:
    |m_peak - m_base| = -2.5 log10 ((A(u0) * f_S + f_L + f_N)/(f_S + f_L + f_N))

    Parameters
    ----------
    u0 : float or array
        Dimensionless source-lens angular separation, closest approach
    f_S : float or array
        Flux from source (arbitrary units)
    f_L : float or array
        Flux from lens (arbitrary units)
    f_N : float or array
        Flux from neighbors (arbitrary units)

    Return
    ------
    m_bump : float or array
        Bump magnitude
    """
    A = calc_magnification(u0)
    m_bump = 2.5 * np.log10((A * f_S + f_L + f_N) / (f_S + f_L + f_N))

    return m_bump


def calc_app_mag(r, M, E, f):
    """
    Calculate the apparent magnitude (i.e. distance modulus
    plus extinction)

    Parameters
    ----------
    M : float or array
        Absolute magnitude of star

    r : float or array
        Distance of star from sun (in kpc)

    E : float or array
        Extinction law

    f : float or array
        Coefficient for that particular band or whatever

    Return
    ------
    m : float or array
        Apparent magnitude of star

    """
    m = calc_DM(r, M) + calc_ext(E, f)

    return m


def calc_DM(r, M):
    """
    Calculate the distance modulus: m = M + 5log10(100*r/kpc)

    Parameters
    ----------
    M : float or array
        Absolute magnitude of star

    r : float or array
        Distance of star from sun (in kpc)

    Return
    ------
    m : float or array
        Apparent magnitude of star

    """
    m = M + 5 * np.log10(100 * r)

    return m


def calc_ext(E, f):
    """
    Calculate the magnitude of extinction.

    Parameters
    ----------
    E : float or array
        Extinction law

    f : float or array
        Coefficient for that particular band or whatever

    Return
    ------
    m_E : float or array
        Magnitude of extinction

    """
    m_E = E * f

    return m_E


def get_Alambda_AKs(red_law_name, lambda_eff):
    """
    Get Alambda/AKs. NOTE: this doesn't work for every law in PopStar!
    Naming convention is not consistent. Change PopStar or add if statements?

    Parameters
    ----------
    red_law_name : str
        The name of the reddening law
    lambda_eff : float
        Wavelength in microns

    Return
    ------
    Alambda_AKs : float
        Alambda/AKs

    """
    red_law_class = getattr(reddening, 'RedLaw' + red_law_name)
    red_law = red_law_class()
    red_law_method = getattr(red_law, red_law_name)
    Alambda_AKs = red_law_method(lambda_eff, 1)

    return Alambda_AKs


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


