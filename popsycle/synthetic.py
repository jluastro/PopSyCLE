#! /usr/bin/env python
"""
synthetic.py
Functions (and their associated functions) for running the PopSyCLE pipeline.
Including:
- write_galaxia_params
- perform_pop_syn
- calc_events
- refine_events
- refine_binary_events
"""
import numpy as np
import h5py
import math
from astropy import units
from scipy.stats import maxwell
import astropy.coordinates as coord
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import Angle  # Angles
from astropy.table import Table, Column, MaskedColumn
from astropy.table import vstack
from spisea.imf import imf
from spisea import synthetic, evolution, reddening, ifmr
from spisea.imf.multiplicity import MultiplicityResolvedDK
from scipy.spatial import cKDTree
import time
import datetime
import gc
import subprocess
import os
from sklearn import neighbors
import itertools
from multiprocessing import Pool, Value, Lock, Manager
import inspect
import numpy.lib.recfunctions as rfn
import copy
from distutils import spawn
from popsycle import ebf
from popsycle.filters import transform_ubv_to_ztf
from popsycle import utils
import astropy.units as unit
from popsycle import orbits
import pandas as pd
from bagle import model
from scipy.signal import find_peaks
from collections import Counter
from operator import itemgetter
from popsycle import binary_utils
from astropy.io import fits

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
# Dictionary for SPISEA IFMR objects
##########
IFMR_dict = {}
IFMR_dict['Raithel18'] = ifmr.IFMR_Raithel18()
IFMR_dict['Spera15'] = ifmr.IFMR_Spera15()
IFMR_dict['SukhboldN20'] = ifmr.IFMR_N20_Sukhbold()
##########

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
photometric_system_dict['ztf'] = ['g', 'r', 'i']

##########
# List of all supported photometric systems and filters with SPISEA labels
##########
all_filt_list = ['ubv,U', 'ubv,B', 'ubv,V', 'ubv,I', 'ubv,R',
                 'ukirt,H', 'ukirt,K', 'ukirt,J']

##########
# List of all supported multiplicity classes
##########
multiplicity_list = {'None': None,
                     'ResolvedDK': MultiplicityResolvedDK}

###########################################################################
############# Population synthesis and associated functions ###############
###########################################################################


def write_galaxia_params(output_root,
                         longitude, latitude, area,
                         seed=None,
                         output_dir='./',
                         photo_sys='UBV',
                         mag_color_names='V,B-V',
                         app_mag_limits=[-1000, 1000],
                         abs_mag_limits=[-1000, 1000],
                         color_limits=[-1000, 1000],
                         geometry_option=1,
                         f_sample=1,
                         pop_ID=-1,
                         warp_flare_on=1,
                         r_max=30):
    """
    Given an object root, sky location and area, creates the parameter
    file that Galaxia requires for running. User can also specify a seed for
    Galaxia to use in its object generation.

    Parameters
    ----------
    output_root : str
        The thing you want the output files to be named
        Examples include
        ``myout``,
        ``/some/path/to/myout``,
        ``../back/to/some/path/myout``.

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    seed : int, optional
         Seed Galaxia will use to generate objects. If not set, script will
         generate a seed from the current time. Setting this seed guarantees
         identical results.

    Returns
    -------
    <output_root>_galaxia_params.txt : str
        A text file with the parameters that Galaxia requires to run.
    """

    if seed is None:
        # Grab last two digits of the current Epochal time
        seed = int(str(time.time())[-2:])

    params = [
        "outputFile %s" % output_root,
        "outputDir %s" % output_dir,
        "photoSys %s" % photo_sys,
        "magcolorNames %s" % mag_color_names,
        "appMagLimits[0] %s" % app_mag_limits[0],
        "appMagLimits[1] %s" % app_mag_limits[1],
        "absMagLimits[0] %s" % abs_mag_limits[0],
        "absMagLimits[1] %s" % abs_mag_limits[1],
        "colorLimits[0] %s" % color_limits[0],
        "colorLimits[1] %s" % color_limits[1],
        "geometryOption %s" % geometry_option,
        "longitude %f" % longitude,
        "latitude %f" % latitude,
        "surveyArea %.5f" % area,
        "fSample %s" % f_sample,
        "popID %s" % pop_ID,
        "warpFlareOn %s" % warp_flare_on,
        "seed %i" % seed,
        "r_max %s" % r_max,
        "starType 0",
        "photoError 0"
    ]

    galaxia_param_fname = '%s_galaxia_params.txt' % output_root

    print('** Generating %s **' % galaxia_param_fname)

    with open(galaxia_param_fname, 'w') as f:
        for param in params:
            f.write(param + '\n')
            print('-- %s' % param)


def _check_run_galaxia(output_root, longitude, latitude, area,
                       galaxia_galaxy_model_filename, seed):
    """
    Check that the inputs to run_galaxia are valid

    Parameters
    -----------
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

    galaxia_galaxy_model_filename : str
        Name of the galaxia galaxy model, as outlined at https://github.com/jluastro/galaxia

    seed : int
         Seed Galaxia will use to generate objects. If not set, script will
         generate a seed from the current time. Setting this seed guarantees
         identical results.
    """

    if not isinstance(output_root, str):
        raise Exception('output_root (%s) must be a string.' % str(output_root))

    if not isinstance(longitude, int):
        if not isinstance(longitude, float):
            raise Exception('longitude (%s) must be an integer or a float.' % str(longitude))

    if not isinstance(latitude, int):
        if not isinstance(latitude, float):
            raise Exception('latitude (%s) must be an integer or a float.' % str(latitude))

    if not isinstance(area, int):
        if not isinstance(area, float):
            raise Exception('area (%s) must be an integer or a float.' % str(area))

    # Check that galaxia is in the executable PATH
    if spawn.find_executable('galaxia') is None:
        raise Exception('galaxia is not an executable currently in $PATH')

    # Confrim that galaxia is the PopSyCLE compliant version
    stdout, _ = utils.execute('galaxia --version')
    galaxia_version = stdout.replace('\n', '').split()[1]
    if galaxia_version != '0.7.2.1':
        raise Exception('galaxia must be version 0.7.2.1 installed from https://github.com/jluastro/galaxia')

    # Check the galaxia_galaxy_model_filename for correct type and existence
    if not isinstance(galaxia_galaxy_model_filename, str):
        raise Exception('galaxia_galaxy_model_filename (%s) must be a string.' % str(galaxia_galaxy_model_filename))

    if not os.path.exists(galaxia_galaxy_model_filename):
        raise Exception('galaxia_galaxy_model_filename (%s) does not exist' % galaxia_galaxy_model_filename)

    # Check that GalaxiaData is a valid galaxia folder
    GalaxiaData = None
    for line in open(galaxia_galaxy_model_filename, 'r'):
        if 'GalaxiaData' in line:
            GalaxiaData = line.replace('\n', '').split()[1]
            break

    if GalaxiaData is None:
        raise Exception('GalaxiaData missing from '
                        'galaxia_galaxy_model_filename (%s)' % galaxia_galaxy_model_filename)

    if not os.path.exists(GalaxiaData):
        raise Exception('GalaxiaData (%s) in galaxia_galaxy_model_filename '
                        '(%s) does not exist' % (GalaxiaData, galaxia_galaxy_model_filename))

    if seed is not None:
        if not isinstance(seed, int):
            raise Exception('seed (%s) must be None or an integer.' % str(seed))


def run_galaxia(output_root, longitude, latitude, area,
                galaxia_galaxy_model_filename,
                seed=None,
                output_dir='./',
                photo_sys='UBV',
                mag_color_names='V,B-V',
                app_mag_limits=[-1000, 1000],
                abs_mag_limits=[-1000, 1000],
                color_limits=[-1000, 1000],
                geometry_option=1,
                f_sample=1,
                pop_ID=-1,
                warp_flare_on=1,
                r_max=30):
    """
    Given an object root, sky location and area, creates the parameter
    file that Galaxia requires for running and executes Galaxia.
    User can also specify a seed for Galaxia to use in its object generation.

    Parameters
    ----------
    output_root : str
        The thing you want the output files to be named
        Examples include
        ``myout``,
        ``/some/path/to/myout``,
        ``../back/to/some/path/myout``.

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    galaxia_galaxy_model_filename : str
        Name of the galaxia galaxy model parameter file,
        as outlined at https://github.com/jluastro/galaxia

    seed : int, optional
         Seed Galaxia will use to generate objects. If not set, script will
         generate a seed from the current time. Setting this seed guarantees
         identical results.
    """

    # Error handling/complaining if input types are not right.
    _check_run_galaxia(output_root, longitude, latitude, area,
                       galaxia_galaxy_model_filename, seed)

    # Writes out galaxia params to disk
    write_galaxia_params(output_root=output_root,
                         longitude=longitude,
                         latitude=latitude,
                         area=area,
                         seed=seed,
                         output_dir=output_dir,
                         photo_sys=photo_sys,
                         mag_color_names=mag_color_names,
                         app_mag_limits=app_mag_limits,
                         abs_mag_limits=abs_mag_limits,
                         color_limits=color_limits,
                         geometry_option=geometry_option,
                         f_sample=f_sample,
                         pop_ID=pop_ID,
                         warp_flare_on=warp_flare_on,
                         r_max=r_max)

    # Execute Galaxia
    cmd = 'galaxia -r %s_galaxia_params.txt %s' % (output_root, galaxia_galaxy_model_filename)
    print('** Executing Galaxia with %s_galaxia_params.txt '
          'and %s **' % (output_root, galaxia_galaxy_model_filename))
    t0 = time.time()
    stdout, stderr = utils.execute(cmd)
    print('** STDOUT **')
    print(stdout)
    print('** STDERR **')
    print(stderr)
    t1 = time.time()
    print('Galaxia complete')
    print('galaxia runtime : {0:f} s'.format(t1 - t0))

    ##########
    # Make log file
    ##########

    stdout, _ = utils.execute('which galaxia')
    galaxia_path = stdout.replace('\n', '')

    now = datetime.datetime.now()
    popsycle_path = os.path.dirname(inspect.getfile(run_galaxia))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=popsycle_path).decode('ascii').strip()
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'longitude , ' + str(longitude) + '\n'
    line2 = 'latitude , ' + str(latitude) + '\n'
    line3 = 'area , ' + str(area) + '\n'
    line3b = 'galaxia_galaxy_model_filename , ' + galaxia_galaxy_model_filename + '\n'
    line4 = 'seed , ' + str(seed) + '\n'

    line8 = 'VERSION INFORMATION' + '\n'
    line9 = str(now) + ' : creation date' + '\n'
    line10 = popsycle_hash + ' : PopSyCLE commit' + '\n'
    line11 = galaxia_path + ' : galaxia path' + '\n'

    line12 = 'OTHER INFORMATION' + '\n'
    line13 = str(t1 - t0) + ' : total runtime (s)' + '\n'

    line17 = 'FILES CREATED' + '\n'
    line18 = output_root + '.ebf : ebf file' + '\n'

    with open(output_root + '_galaxia.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, line3b, line4,
                        empty_line, line8, dash_line, line9, line10, line11,
                        empty_line, line12, dash_line, line13,
                        empty_line, line17, dash_line, line18])


def _check_perform_pop_syn(ebf_file, output_root, iso_dir,
                           IFMR,
                           bin_edges_number,
                           BH_kick_speed_mean, NS_kick_speed_mean,
                           additional_photometric_systems,
                           multiplicity, binning,
                           overwrite, seed,
                           n_proc, verbose):
    """
    Checks that the inputs of perform_pop_syn are valid

    Parameters
    -----------
    ebf_file : str or ebf file
        str : name of the ebf file from Galaxia

    output_root : str
        The thing you want the output files to be named
        Examples:
           'myout'
           '/some/path/to/myout'
           '../back/to/some/path/myout'

    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)

    IFMR : string
        The name of the IFMR object from SPISEA. For more information on these objects see ifmr.py
        in SPISEA. 
        'Raithel18' = IFMR_Raithel18
        'Spera15' = IFMR_Spera15
        'SukhboldN20' = IFMR_N20_Sukhbold

    bin_edges_number : int
        Number of edges for the bins
            bins = bin_edges_number - 1
        Total number of bins is
            N_bins = (bin_edges_number - 1)**2
        If set to None (default), then number of bins is
            bin_edges_number = int(60 * 2 * radius) + 1

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
        Galaxia / SPISEA's ubv photometry and appended to the output files.

    multiplicity: object or None
        If a resovled multiplicity object is specified, 
        the table will be generated with resolved multiples.
    
    binning : bool
        If set to True, bins files as specified by bin_edges_numbers or default.
        If set to False, no bins (SET TO FALSE IF DOING FULL SKY DOWNSAMPLED).
    
    overwrite : bool
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.

    seed : int or None
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for SPISEA and PopSyCLE.
        
    n_proc : int
        Number of processors to use in the parallel processing. Note that
        calculations are memory intensive and n_proc > few should likely not be
        used to avoid memory overruns that may result in data corruption.

    verbose : int
        Level of debugging information to print to stderr. Set to 0 for minimal
        information. Coarse timing at 2 and fine timing at 4.
    """

    if not isinstance(ebf_file, str):
        raise Exception('ebf_file (%s) must be a string.' % str(ebf_file))

    if ebf_file[-4:] != '.ebf':
        raise Exception('ebf_file (%s) must be an ebf file.' % str(ebf_file))

    if not isinstance(output_root, str):
        raise Exception('output_root (%s) must be a string.' % str(output_root))

    if not isinstance(iso_dir, str):
        raise Exception('iso_dir (%s) must be a string.' % str(iso_dir))

    if not os.path.exists(iso_dir):
        raise Exception('iso_dir (%s) must exist' % str(iso_dir))

    if not isinstance(IFMR, str):
        raise Exception('IFMR (%s) must be a string.' % str(IFMR))

    if IFMR not in IFMR_dict:
        exception_str = '(%s) is not a valid IFMR. ' \
                        'IFMR must be a key in ' \
                        'IFMR_dict. \n' \
                        'Acceptable values are the strings : ' %str(IFMR)
        for IFMR in IFMR_dict:
            exception_str += '%s, ' % IFMR
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    if bin_edges_number is not None:
        if not isinstance(bin_edges_number, int):
            raise Exception('bin_edges_number (%s) must be None or an integer.' % str(bin_edges_number))

    if not isinstance(BH_kick_speed_mean, int):
        if not isinstance(BH_kick_speed_mean, float):
            raise Exception('BH_kick_speed_mean (%s) must be an integer or a float.' % str(BH_kick_speed_mean))

    if not isinstance(NS_kick_speed_mean, int):
        if not isinstance(NS_kick_speed_mean, float):
            raise Exception('NS_kick_speed_mean (%s) must be an integer or a float.' % str(NS_kick_speed_mean))

    if multiplicity is not None:
        if not isinstance(multiplicity, MultiplicityResolvedDK):
            raise Exception('multiplicity must be None or a subclass of MultiplicityResolvedDK.')
    
    if not isinstance(binning, bool):
        raise Exception('binning (%s) must be a boolean.' % str(binning))
    
    if not isinstance(overwrite, bool):
        raise Exception('overwrite (%s) must be a boolean.' % str(overwrite))

    if seed is not None:
        if not isinstance(seed, int):
            raise Exception('seed (%s) must be None or an integer.' % str(seed))
    
    if not isinstance(n_proc, int):
        raise Exception('n_proc (%s) must be an integer.' % str(seed))
    
    if not isinstance(verbose, int):
        raise Exception('verbose (%s) must be an integer.' % str(seed))
    

    if additional_photometric_systems is not None:
        if not isinstance(additional_photometric_systems, list):
            raise Exception('additional_photometric_systems (%s) must '
                            'either None or a list of strings.' % str(additional_photometric_systems))

        for photometric_system in additional_photometric_systems:
            if photometric_system not in photometric_system_dict:
                exception_str = 'photometric_system must be a key in ' \
                                'photometric_system_dict. \n' \
                                'Acceptable values are : '
                for photometric_system in photometric_system_dict:
                    exception_str += '%s, ' % photometric_system
                exception_str = exception_str[:-2]
                raise Exception(exception_str)


def perform_pop_syn(ebf_file, output_root, iso_dir,
                    IFMR,
                    bin_edges_number=None,
                    BH_kick_speed_mean=50, NS_kick_speed_mean=400,
                    additional_photometric_systems=None,
                    multiplicity=None, binning = True,
                    overwrite=False, seed=None,
                    n_proc=1, verbose=0):
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
        Examples include 
        ``myout``, 
        ``/some/path/to/myout``
        ``../back/to/some/path/myout``.

    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)

    IFMR : string
        The name of the IFMR object from SPISEA. For additional information on these objects see ifmr.py
        in SPISEA. 
        'Raithel18' = IFMR_Raithel18
        'Spera15' = IFMR_Spera15
        'SukhboldN20' = IFMR_N20_Sukhbold

    bin_edges_number : int, optional
        Number of edges for the bins such that
        ``bins = bin_edges_number - 1``.
        The Total number of bins is
        ``N_bins = (bin_edges_number - 1)**2``.
        If set to None (default), then number of bins is
        ``bin_edges_number = int(60 * 2 * radius) + 1``

    BH_kick_speed_mean : float, optional
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.

    NS_kick_speed_mean : float, optional
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract

    additional_photometric_systems : list of strs, optional
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.

    multiplicity: object
        If a resovled multiplicity object is specified, 
        the table will be generated with resolved multiples.
        Default is None.
    
    binning : bool
        If set to True, bins files as specified by bin_edges_numbers or default.
        If set to False, no bins (SET TO FALSE IF DOING FULL SKY DOWNSAMPLED).
        Default is True.

    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.
        Default is False.

    seed : int, optional
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for SPISEA and PopSyCLE.
        Default None.

    n_proc : int
        Number of processors to use in the parallel processing. Note that
        calculations are memory intensive and n_proc > few should likely not be
        used to avoid memory overruns that may result in data corruption.
        Default is 1.

    verbose : int
        Level of debugging information to print to stderr. Set to 0 for minimal
        information. Coarse timing at 2 and fine timing at 4.
        Default is 0.

    Returns
    -------
    <output_root>.h5 : hdf5 file
        NOTE: This is what _bin_lb_hdf5 returns.
        A compund type hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects
        and stars sorted into those bins.

    <output_root>_label.fits : fits file
        NOTE: This is what make_label_file returns.
        A fits file that shows the correspondence between dataset name,
        latitude, longitude, and number of objects in that bin.
    """
    # Check whether files exist
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
    _check_perform_pop_syn(ebf_file, output_root, iso_dir,
                           IFMR,
                           bin_edges_number,
                           BH_kick_speed_mean, NS_kick_speed_mean,
                           additional_photometric_systems,
                           multiplicity, binning,
                           overwrite, seed,
                           n_proc, verbose)

    ##########
    # Start of code
    #########

    # Set random seed
    np.random.seed(seed)

    t0 = time.time()

    #########
    # Read in only what you need
    # Control yourself... take only what you need from it
    # i.e. the log file, popid, metallicity, and age
    ########
    t = ebf.read_ind(ebf_file, '/log', 0)
    popid_array = ebf.read(ebf_file, '/popid')
    metallicity_array = ebf.read(ebf_file, '/feh')
    age_array = ebf.read(ebf_file, '/age')

    # Just some counting/error checking things
    n_stars = len(popid_array)  # How many stars from Galaxia
    co_counter = 0  # Number of compact objects made
    unmade_cluster_counter = 0 # Number of clusters that are too small and get skipped over
    unmade_cluster_mass = 0 # Current mass in the unmade clusters

    # Convert log to useful dictionary.
    ebf_log = make_ebf_log(t)

    ##########
    #  Fetch the center point and area of the survey.
    ##########
    b = float(ebf_log['latitude'])
    l = float(ebf_log['longitude'])
    surveyArea = float(ebf_log['surveyArea'])

    ##########
    # Make the bins for the latitude and longitude. All the information
    # will be sorted into each of these bins in order to
    # handle large arrays, etc.
    ##########
    bin_edges_number, lat_bin_edges, long_bin_edges = _get_bin_edges(l, b,
                                                                     surveyArea,
                                                                     bin_edges_number)

    ##########
    # Create h5py file to store lat/long binned output
    ##########
    h5file = h5py.File(output_root + '.h5', 'w')
    if binning == True:
        h5file['lat_bin_edges'] = lat_bin_edges
        h5file['long_bin_edges'] = long_bin_edges
    if 'galaxyModelFile' in ebf_log:
        h5file['galaxyModelFile'] = ebf_log['galaxyModelFile']
    h5file.close()
    
    # Make one for companions if multiplicity
    if multiplicity is not None:
        h5file_comp = h5py.File(output_root + '_companions' + '.h5', 'w')
        if binning:
            h5file_comp['lat_bin_edges'] = lat_bin_edges
            h5file_comp['long_bin_edges'] = long_bin_edges
        if 'galaxyModelFile' in ebf_log:
            h5file_comp['galaxyModelFile'] = ebf_log['galaxyModelFile']
        h5file_comp.close()
        
    ##########
    # Reassign ages for stars that are less than logage 5.01 
    # or greater than logage 10.14, since those are the limits of
    # SPISEA evolution. Justification: in writeup/paper.
    ##########
    young_id = np.where(age_array <= 5.01)[0]
    age_array[young_id] = 5.0101
    old_id = np.where(age_array >= 10.14)[0]
    age_array[old_id] = 10.1399

    # Initialize a running counter for the "unique ID" for both
    # stars and compact objects. Note these must be globals to work
    # across processes.
    next_id_stars_val_in = Value('i', 0)
    next_id_co_val_in = Value('i', n_stars)

    # Single-Threaded
    # _mp_init_worker(mp_lock, next_id_stars_val_in, next_id_co_val_in)

    # Multi-Threaded
    # Setup up a pool for multiprocessing.
    mp_lock = Lock()
    mp_pool = Pool(processes=n_proc,
                   initializer=_mp_init_worker,
                   initargs=(mp_lock, next_id_stars_val_in, next_id_co_val_in))

    ##########
    # Loop through population ID (i.e. bulge, disk, halo, ...) in order to
    # design bin sizes appropriate for each population. Dividing by population
    # is convenient because each has very different
    # age ranges and radius ranges.
    ##########
    for pid in range(10):
        print('******************** Starting popid ' + str(pid))
        if np.sum(popid_array == pid) == 0:
            print('No stars with this pid. Skipping!')
            continue

        popid_idx = np.where(popid_array == pid)[0]

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
            print(f'********** Starting log age bin {logt_bins[aa]:.1f}')
            # Mid-point age of bin.
            age_of_bin = (logt_bins[aa] + logt_bins[aa + 1]) / 2.0

            # Fetch the stars in this age bin.
            age_idx = np.where((age_array[popid_idx] >= logt_bins[aa]) &
                               (age_array[popid_idx] < logt_bins[aa + 1]))[0]

            # Create a KDTree from randomly selected stars in the
            # pop_id / age_bin used for calculating extinction to luminous
            # white dwarfs. Because the same KDTree is used for each age-bin,
            # two compact objects randomly selected to have nearly identical
            # positions would have identical extinctions. This low
            # probability event is a reasonable trade-off for the reduced
            # compute time gained by only constructing the KDTree once.
            # The size of this KDtree is set by the need to sample the
            # spatial volume densly enough. We estimate ~2 million stars per box
            # should suffice; however, whole-sky runs might need more and currently
            # may have less accurate extinctions.
            exbv_arr4kdt, kdt_star_p = _make_extinction_kdtree(ebf_file, popid_idx[age_idx])

            if IFMR == 'Raithel18':
                # Only run at solar metallicity for Raithel18
                # -99 to 99 will capture all possible values of metallicity, but the Galaxia range is closer to -2 to 2
                feh_bins = [-99, 99]
                feh_vals = [0.0]
                
            elif IFMR == 'Spera15' or 'SukhboldN20':
                # Break into 4 hardcoded metallicity bins for Spera15 and SukhboldN20
                # By starting at -99 and ending at 99 we will capture all possible metallicty values, but the true distribution is much narrower
                feh_bins = [-99, -1.279, -0.500, 0.00, 99]
                feh_vals = [-1.39, -0.89, -0.25, 0.30]

            mp_proc = [] # This will hold the multiprocessing task results for this age bin.

            ##########
            # Loop through metallicity bins
            ##########
            for bb in range(len(feh_bins) - 1):
                print('***** Starting metallicity bin ', feh_vals[bb])
                # Value of the metallicity bin
                metallicity_of_bin = feh_vals[bb]

                # Fetch the stars in this metallicity bin
                feh_idx = np.where((metallicity_array[popid_idx[age_idx]] >= feh_bins[bb]) &
                                   (metallicity_array[popid_idx[age_idx]] < feh_bins[bb + 1]))[0]
                len_adx = len(feh_idx)

                # Figure out how many bins we will need.
                #   -- breaking up into managable chunks of 50,000 stars each.
                #   -- https://docs.google.com/document/d/1Q2i5I7s-tZYH0lsWOEb1RH6_sv76cDm7TtPN6c7RSII/edit#
                #      for detailed notes on this choice.
                num_stars_in_bin = 5e4
                num_bins = int(math.ceil(len_adx / num_stars_in_bin))

                ##########
                # Loop through bins of 50,000 stars at a time.
                ##########
                for nn in range(num_bins):
                    n_start = int(nn * num_stars_in_bin)
                    n_stop = int((nn + 1) * num_stars_in_bin)

                    bin_idx = popid_idx[age_idx[feh_idx[n_start:n_stop]]]
                    print(f'** Starting sub-bin {nn:2d} with {len(bin_idx)} stars in bin')

                    # # Single-threaded version. - START
                    # count_out = _process_popsyn_stars_in_bin(bin_idx, age_of_bin, metallicity_of_bin,
                    #                                          popid_array, age_array, lat_bin_edges, long_bin_edges,
                    #                                          ebf_file,
                    #                                          kdt_star_p, exbv_arr4kdt,
                    #                                          iso_dir, IFMR, NS_kick_speed_mean, BH_kick_speed_mean,
                    #                                          multiplicity,
                    #                                          additional_photometric_systems,
                    #                                          t0, binning, seed, output_root, verbose=verbose)
                    #
                    # co_count_new, unmade_cluster_count_new, unmade_cluster_mass_new = count_out
                    #
                    # # Update running counters.
                    # co_counter += co_count_new
                    # unmade_cluster_counter += unmade_cluster_count_new
                    # unmade_cluster_mass += unmade_cluster_mass_new
                    # # Single-threaded version. - STOP

                    # Multi-threaded version. - START
                    args = (bin_idx, age_of_bin, metallicity_of_bin,
                            popid_array, age_array, lat_bin_edges, long_bin_edges,
                            ebf_file,
                            kdt_star_p, exbv_arr4kdt,
                            iso_dir, IFMR, NS_kick_speed_mean, BH_kick_speed_mean,
                            multiplicity,
                            additional_photometric_systems,
                            t0, binning, seed, output_root, verbose)

                    # proc_tmp = Process(target=_process_popsyn_stars_in_bin, args=args)
                    mp_res_tmp = mp_pool.apply_async(_process_popsyn_stars_in_bin, args=args)

                    mp_proc.append(mp_res_tmp)

            # Complete this age bin.
            # Gather results from parallel runs over this whole age_bin
            results = [proc.get() for proc in mp_proc]
            
            if len(results) > 0:
                co_counter += sum(list(zip(*results))[0])
                unmade_cluster_counter += sum(list(zip(*results))[1])
                unmade_cluster_mass += sum(list(zip(*results))[2])
            # Multi-threaded version. - STOP


            del kdt_star_p, exbv_arr4kdt
            gc.collect()

    # Multi-Threaded: clean up pool
    mp_pool.close()
    mp_pool.join()
    
    t1 = time.time()
    print('perform_pop_syn runtime : {0:f} s'.format(t1 - t0))

    ##########
    # Figure out how much stuff got binned.
    ##########
    binned_counter = 0
    hf = h5py.File(output_root + '.h5', 'r')
    for key in hf:
        if 'bin_edges' not in key and 'galaxyModel' not in key:
            binned_counter += len(hf[key])

    ##########
    # Make label file containing information about l,b bins
    ##########
    if binning == True:
        make_label_file(output_root, overwrite=overwrite)
    else:
        make_label_file_no_bins(output_root, overwrite=overwrite)

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    popsycle_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=popsycle_path).decode('ascii').strip()
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
    line7 = 'IFMR , ' + str(IFMR) + '\n'
    line8 = 'seed , ' + str(seed) + '\n'

    line9 = 'VERSION INFORMATION' + '\n'
    line10 = str(now) + ' : creation date' + '\n'
    line11 = popstar_hash + ' : SPISEA commit' + '\n'
    line12 = popsycle_hash + ' : PopSyCLE commit' + '\n'

    line13 = 'OTHER INFORMATION' + '\n'
    line14 = str(t1 - t0) + ' : total runtime (s)' + '\n'
    line15 = str(n_stars) + ' : total stars from Galaxia' + '\n'
    line16 = str(co_counter) + ' : total compact objects made' + '\n'
    line17 = str(binned_counter) + ' : total things binned' + '\n'
    line18 = str(unmade_cluster_counter) + ' : # of unmade clusters (< 100 M_sun)' + '\n'
    line19 = str(unmade_cluster_mass) + ' : # mass in unmade clusters (M_sun)' + '\n'

    line20 = 'FILES CREATED' + '\n'
    line21 = output_root + '.h5 : HDF5 file' + '\n'
    line22 = output_root + '_label.fits : label file' + '\n'
    if multiplicity is not None:
        line23 = output_root + '_companions.h5 : HDF5 file' + '\n'
    else:
        line23 = 'No companions h5 file as multiplicity = None' + '\n'

    with open(output_root + '_perform_pop_syn.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, line4, line5,
                        line6, line7, line8, empty_line, line9, dash_line,
                        line10, line11, line12, empty_line,line13, dash_line,
                        line14, line15, line16, line17, line18, line19, empty_line, line20, dash_line,
                        line21, line22, line23])

    ##########
    # Informative print statements.
    ##########
    if (n_stars + co_counter) != binned_counter:
        print('***************** WARNING ******************')
        print('Number of things in != number of things out.')
        print('********************************************')

    print('******************** INFO **********************')
    print('Total number of stars from Galaxia: ' + str(n_stars))
    print('Total number of compact objects made: ' + str(co_counter))
    print('Total number of things binned: ' + str(binned_counter))

    return

def _mp_init_worker(lock_in, next_id_stars_val_in, next_id_co_val_in):
    """
    Initialize the multiprocessing pool and save global tracking variables.
    Saves global variables: lock, next_id_stars_val, next_id_co_val

    Parameters
    -----------
    lock_in : multiprocessing.Lock object
    next_id_stars_val_in : multiprocessing.Value object
    next_id_co_val_in : multiprocessing.Value object

    """
    global lock, next_id_stars_val, next_id_co_val

    lock = lock_in
    next_id_stars_val = next_id_stars_val_in
    next_id_co_val = next_id_co_val_in

    return


def _process_popsyn_stars_in_bin(bin_idx, age_of_bin, metallicity_of_bin,
                                 popid_array, age_array, lat_bin_edges, long_bin_edges,
                                 ebf_file, kdt_star_p, exbv_arr4kdt,
                                 iso_dir, IFMR, NS_kick_speed_mean, BH_kick_speed_mean,
                                 multiplicity, additional_photometric_systems, t0,
                                 binning, seed,
                                 output_root, verbose=0):
    
    """
    Processess the objects in a bin including 
    creating the star dict, generating the compact objects
    and companions (if applicable) and writing to hdf5 file.
    
    Parameters
    -----------
    bin_idx : array
        Array of indices of stars that are in the bin
    
    age_of_bin : float
        Mid-point age of bin
        
    metallicity_of_bin : float
        Value of the metallicity bin
        
    pop_id_array : array
        Array of population ids (i.e. 0-9) of Galaxia stars
    
    age_array : array
        Array of ages of Galaxia stars
    
    lat_bin_edges : array
        Edges for the latitude binning (deg)

    long_bin_edges : array
        Edges for the longitude binning (deg)
    
    ebf_file : str or ebf file
        str : name of the ebf file from Galaxia
        ebf file : actually the ebf file from Galaxia
        
    kdt_star_p : scipy KDTree
        KDTree of the momenta of a random selection of Galaxia stars.
        
    exbv_arr4kdt : array
        Exteniction of the stars in the KD Tree.
    
    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)
        
    IFMR : string
        The name of the IFMR object from SPISEA. For additional information on these objects see ifmr.py
        in SPISEA. 
        'Raithel18' = IFMR_Raithel18
        'Spera15' = IFMR_Spera15
        'SukhboldN20' = IFMR_N20_Sukhbold
        
    NS_kick_speed_mean : float
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract
        
    BH_kick_speed_mean : float
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.
        
    multiplicity : multiplicity: object
        If a resovled multiplicity object is specified, 
        the table will be generated with resolved multiples.
        Default is None.
        
    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.
    
    t0 : float
        Start time for timing purposes
    
    binning : bool
        If set to True, bins files as specified by bin_edges_numbers or default.
        If set to False, no bins (SET TO FALSE IF DOING FULL SKY DOWNSAMPLED).
        Default is True.

    seed : int, optional
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for SPISEA and PopSyCLE.
        Default None.
        
    output_root :str
        The thing you want the output files to be named
        Examples:
           'myout'
           '/some/path/to/myout'
           '../back/to/some/path/myout'
    
    verbose : int, optional
        Level of debugging information to print to stderr. Set to 0 for minimal
        information. Coarse timing at 2 and fine timing at 4.
        Default is 0.
    
    Returns
    --------
    n_co_new : int
        Number of new compact objects generated per bin.
    
    unmade_cluster_counter_new : int
        Updated number of unmade clusters (<= 100 M_sun)

    unmade_cluster_mass_new: float
        The current mass in the unmade clusters (<= 100 M_sun)
    
    """
    # Load up our global (thread-safe) variables.
    global lock, next_id_stars_val, next_id_co_val

    ##########
    # Fill up star_dict
    ##########
    star_dict = {}

    # Multi-Threaded (works in Single-Threaded too)
    with lock:
        star_dict['obj_id'] = np.arange(len(bin_idx)) + next_id_stars_val.value
        next_id_stars_val.value += len(bin_idx)


    star_dict['age'] = age_array[bin_idx]
    star_dict['popid'] = popid_array[bin_idx]

    _load_galaxia_into_star_dict(star_dict,
                                 bin_idx,
                                 ebf_file,
                                 additional_photometric_systems)
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
    mass_in_bin = np.sum(star_dict['mass'])
    stars_in_bin = {}
    for key, val in star_dict.items():
        stars_in_bin[key] = val
    cluster_tmp, unmade_cluster_counter_new, unmade_cluster_mass_new = _make_cluster(iso_dir=iso_dir,
                                                                                     log_age=age_of_bin,
                                                                                     currentClusterMass=mass_in_bin,
                                                                                     multiplicity=multiplicity,
                                                                                     IFMR=IFMR,
                                                                                     feh=metallicity_of_bin, seed=seed,
                                                                                     additional_photometric_systems=additional_photometric_systems)
    # Make a temporary table of the compact object primaries. Will eventually be appended to star_dict.
    co_dict, next_id_co = _make_co_dict(age_of_bin,
                                        cluster_tmp,
                                        stars_in_bin,
                                        kdt_star_p, exbv_arr4kdt,
                                        BH_kick_speed_mean=BH_kick_speed_mean,
                                        NS_kick_speed_mean=NS_kick_speed_mean,
                                        additional_photometric_systems=additional_photometric_systems,
                                        multiplicity=multiplicity,
                                        seed=seed)

    #########
    # If there are multiples make companions table
    #########
    if multiplicity is not None:
        star_dict, companions_table = _make_companions_table(cluster=cluster_tmp,
                                                             star_dict=star_dict,
                                                             co_dict=co_dict,
                                                             additional_photometric_systems=additional_photometric_systems,
                                                             t0=t0, verbose=verbose)

        # Save companion table
        with lock:
            if companions_table is not None:
                # Save and bin in l, b all companions
                # Or just save if binning = False
                if binning == True:
                    if co_dict is not None:
                        _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                                     co_dict, output_root + '_companions',
                                     companion_obj_arr=companions_table)
                    _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                                 stars_in_bin, output_root + '_companions',
                                 companion_obj_arr=companions_table)

                else:
                    if co_dict is not None:
                        _no_bins_hdf5(co_dict, output_root + '_companions',
                                      companion_obj_arr=companions_table)
                    _no_bins_hdf5(stars_in_bin, output_root + '_companions',
                                  companion_obj_arr=companions_table)

    if verbose > 3: print(f'test5 {time.time() - t0:.2f} sec')
    del cluster_tmp

    ##########
    #  Save and bin in in l, b all stars and compact objects.
    #  Or just save if binning = False
    ##########
    n_co_new = 0
    with lock:
        if binning == True:
            if co_dict is not None:
                n_co_new = len(co_dict['mass'])
                _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                             co_dict, output_root)
            _bin_lb_hdf5(lat_bin_edges, long_bin_edges,
                         stars_in_bin, output_root)

        else:
            if co_dict is not None:
                n_co_new = len(co_dict['mass'])
                _no_bins_hdf5(co_dict, output_root)
            _no_bins_hdf5(stars_in_bin, output_root)

    if verbose > 1: print(f'test6 {time.time() - t0:.2f} sec')

    ##########
    # Done with galaxia output in dictionary t and ebf_log.
    # Garbage collect in order to save space.
    ##########
    del star_dict


    return n_co_new, unmade_cluster_counter_new, unmade_cluster_mass_new


def _make_extinction_kdtree(ebf_file, indices, max_num_stars_in_kdt_p=int(2e6)):
    """
    Using the stars at the specified indices (but not more than
    2 million of them), construct a KD tree on px, py, pz and also
    note the extinction for each of these same stars. Later,
    we will use this to make a new star at some position,
    query for the nearest star in position in the KD tree and
    assign the associated extinction.

    Parameters
    -----------
    ebf_file : str
        The EBF file from Galaxia that contains all the stars.

    indices : array
        The indices of the stars in this EBF file used to construct
        the KD tree. Note, if the indices array is too large, then
        we will randomly choose a smaller subset of these indices for
        KDTree construction.

    Returns
    -------
    kdt_star_p : KDTree
        The KDTree of Galaxia stars on position [px, py, pz].
    exbv_arr4kdt : numpy.array
        The array of associated extinction values for the stars
        in the KDTree.
    """
    kdt_star_p = None
    exbv_arr4kdt = None

    len_idx = len(indices)

    if len_idx > 0:
        num_kdtree_samples = int(min(len_idx, max_num_stars_in_kdt_p))
        kdt_idx = np.random.choice(np.arange(len_idx),
                                   size=num_kdtree_samples,
                                   replace=False)
        bin_idx = indices[kdt_idx]
        star_px = ebf.read_ind(ebf_file, '/px', bin_idx)
        star_py = ebf.read_ind(ebf_file, '/py', bin_idx)
        star_pz = ebf.read_ind(ebf_file, '/pz', bin_idx)
        star_xyz = np.array([star_px, star_py, star_pz]).T
        kdt_star_p = cKDTree(star_xyz)
        exbv_arr4kdt = ebf.read_ind(ebf_file, '/exbv_schlegel', bin_idx)
        del bin_idx, star_px, star_py, star_pz

    return exbv_arr4kdt, kdt_star_p


def _load_galaxia_into_star_dict(star_dict, bin_idx, ebf_file, additional_photometric_systems):
    """
    The dictionary is edited in place, so nothing is returned.
    
    Parameters
    -----------
    star_dict : dictionary
        The number of entries for each key is the number of stars.
    
    bin_idx : array
        Array of indices of stars that are in the bin
        
    ebf_file : str or ebf file
        str : name of the ebf file from Galaxia
        ebf file : actually the ebf file from Galaxia
        
    additional_photometric_systems : list of strs
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.
    """
    star_dict['zams_mass'] = ebf.read_ind(ebf_file, '/smass', bin_idx)
    star_dict['mass'] = ebf.read_ind(ebf_file, '/mact', bin_idx)
    star_dict['systemMass'] = copy.deepcopy(star_dict['mass'])
    star_dict['px'] = ebf.read_ind(ebf_file, '/px', bin_idx)
    star_dict['py'] = ebf.read_ind(ebf_file, '/py', bin_idx)
    star_dict['pz'] = ebf.read_ind(ebf_file, '/pz', bin_idx)
    star_dict['vx'] = ebf.read_ind(ebf_file, '/vx', bin_idx)
    star_dict['vy'] = ebf.read_ind(ebf_file, '/vy', bin_idx)
    star_dict['vz'] = ebf.read_ind(ebf_file, '/vz', bin_idx)
    star_dict['exbv'] = ebf.read_ind(ebf_file, '/exbv_schlegel', bin_idx)
    star_dict['glat'] = ebf.read_ind(ebf_file, '/glat', bin_idx)
    star_dict['glon'] = ebf.read_ind(ebf_file, '/glon', bin_idx)
    star_dict['mbol'] = ebf.read_ind(ebf_file, '/lum', bin_idx)
    star_dict['grav'] = ebf.read_ind(ebf_file, '/grav', bin_idx)
    star_dict['teff'] = ebf.read_ind(ebf_file, '/teff', bin_idx)
    star_dict['feh'] = ebf.read_ind(ebf_file, '/feh', bin_idx)
    star_dict['rad'] = ebf.read_ind(ebf_file, '/rad', bin_idx)
    star_dict['isMultiple'] = np.zeros(len(bin_idx), dtype=int)
    star_dict['N_companions'] = np.zeros(len(bin_idx), dtype=int)
    star_dict['rem_id'] = np.zeros(len(bin_idx))
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

            del ubv_b, ubv_v, ubv_r, ubv_i, ztf_g, ztf_r, ztf_i


def _get_bin_edges(l, b, surveyArea, bin_edges_number):
    """
    
    Parameters
    -----------    
    l : float
        Galactic longitude, ranging from -180 degrees to 180 degrees
        
    b : float
        Galactic latitude, ranging from -90 degrees to 90 degrees
        
    surveyArea : float
        Area of the sky that will be generated, in square degrees
        
    bin_edges_number :int
        Number of edges for the bins
            bins = bin_edges_number - 1
        Total number of bins is
            N_bins = (bin_edges_number - 1)**2
        If set to None (default in perform_pop_syn), then number of bins is
            bin_edges_number = int(60 * 2 * radius) + 1
    
    Returns
    --------
    bin_edges_number : int
        Modified number of edges for the bins.
        If bin_edges_number = None
            The widths fo the bins are set to 1/2 an arcmin
        Bin numbers are set to at least 3 bins.
        
    lat_bin_edges : array
        Edges for the latitude binning (deg)

    long_bin_edges : array
        Edges for the longitude binning (deg)
    """
    # Extend the edges a bit, that's what the * 1.1 is for
    # (to potentially catch any edge cases.)
    # make bins of size ~1/2 arcmin
    radius = np.sqrt(surveyArea / np.pi)  # degrees
    # Define bin_edges_number, if not given in input.
    if bin_edges_number is None:
        # set the widths to 1/2 arcmin
        bin_edges_number = int(60 * 2 * radius) + 1
    # Make sure we have enough bin edges (min is 3)
    bin_edges_number = max(bin_edges_number, 3)
    # Make sure we have don't have too many bin edges (max is 40)
    bin_edges_number = min(bin_edges_number, 40)
    lat_bin_edges = np.linspace(b - 1.1 * radius, b + 1.1 * radius,
                                bin_edges_number)
    long_bin_edges = np.linspace(l - 1.1 * radius, l + 1.1 * radius,
                                 bin_edges_number)
    # Angle wrapping for longitude
    wrap_id = np.where(long_bin_edges > 180)[0]
    long_bin_edges[wrap_id] -= 360
    return bin_edges_number, lat_bin_edges, long_bin_edges


def _make_co_dict(log_age,
                  cluster,
                  star_dict,
                  kdt_star_p, exbv_arr4kdt,
                  BH_kick_speed_mean=50, NS_kick_speed_mean=400,
                  additional_photometric_systems=None,
                  multiplicity=None,
                  seed=None):
    
    """
    Perform population synthesis.

    Parameters
    -----------
    log_age : float
        log(age/yr) of the cluster you want to make, rounds to nearest 0.01

    cluster : SPISEA.ResolvedCluster
        A SPISEA cluster to use to inject compact objects.

    star_dict : dictionary
        The number of entries for each key is the number of stars.

    kdt_star_p : scipy cKDTree
        KDTree constructed from the positions of randomly selected stars
        that all share the same popid and similar log_age.

    exbv_arr4kdt : numpy array
        Array of galactic extinctions for the stars in kdt_star_p

    BH_kick_speed_mean : float, optional
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.

    NS_kick_speed_mean : float, optional
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract

    additional_photometric_systems : list of strs or None, optional
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.

    seed : int or None, optional
         Seed used to sample the kde tree. If set to any number,
         SPISEA will also be forced to use 42 as a
         random seed for calls to ResolvedCluster.
         Default is None.
         
    multiplicity: object or None, optional
        If a resovled multiplicity object is specified, 
        the table will be generated with resolved multiples.
        Default is None.

    Returns
    -------
    co_dict : dictionary
        Keys are the same as star_dict, just for compact objects.

    next_id_co : int
        Updated next unique ID number (int) that will be assigned to
        the new compact objects created.

    """
    global lock, next_id_co_val

    co_dict = None

    if cluster:
        output = cluster.star_systems
        
        # Create the SPISEA table with just compact objects

        # The compact IDs are:
        # 101: WD, 102: NS, 103: BH
        # -1: PMS, 0: MS, 2: RGB, 3: CHeB, 4: EAGB, 5: TPAGB, 9: WR
        # WE CONSIDER EVERYTHING FROM THE MIST MODELS
        # WITH PHASE 6 (postAGB) TO BE WHITE DWARFS
        compact_ID = np.where((output['phase'] == 101) |
                              (output['phase'] == 102) |
                              (output['phase'] == 103))[0]
        co_table = output[compact_ID]

        # Removes unused columns to conserve memory.
        keep_columns = ['zams_mass', 'isMultiple','systemMass', 'phase', 'mass', 'm_ubv_I', 'm_ubv_R',
                        'm_ubv_B', 'm_ubv_U', 'm_ubv_V', 'm_ukirt_H',
                        'm_ukirt_J', 'm_ukirt_K']
        if multiplicity is not None:
            keep_columns.append('N_companions')
            
        if additional_photometric_systems is not None:
            if 'ztf' in additional_photometric_systems:
                keep_columns += ['m_ztf_g', 'm_ztf_r', 'm_ztf_i']
        co_table.keep_columns(keep_columns)

        # Fill out the rest of co_dict
        if len(co_table['zams_mass']) > 0:

            # Turn astropy table into dictionary to conserve memory.
            co_dict = {}
            co_dict['mass'] = co_table['mass'].data
            co_dict['rem_id'] = co_table['phase'].data
            co_dict['zams_mass'] = co_table['zams_mass'].data
            co_dict['isMultiple'] = co_table['isMultiple'].data
            # makes sure the system mass is current mass
            co_dict['systemMass'] = copy.deepcopy(co_dict['mass'])
            if multiplicity is not None:
                co_dict['N_companions'] = co_table['N_companions'].data

            ##########
            # Assign spherical positions and velocities to all compact objects.
            ##########
            # Create velocities, and positions using Kernel Density Estimation
            # Galactic position (r, l, b)
            # and heliocentric velocity (vx, vy, vz)
            kde_in_data = np.array([star_dict['rad'], star_dict['glat'],
                                    star_dict['glon'], star_dict['vx'],
                                    star_dict['vy'], star_dict['vz']]).T
            # .T because kde.fit needs the rows/columns this way.
            kde = neighbors.KernelDensity(bandwidth=0.0001)
            kde.fit(kde_in_data)

            kde_out_data = kde.sample(len(co_dict['mass']),
                                      random_state=seed)
            co_dict['rad'] = kde_out_data[:, 0]
            co_dict['glat'] = kde_out_data[:, 1]
            co_dict['glon'] = kde_out_data[:, 2]
            co_dict['vx'] = kde_out_data[:, 3]
            co_dict['vy'] = kde_out_data[:, 4]
            co_dict['vz'] = kde_out_data[:, 5]

            ##########
            # Assign x, y, z position.
            ##########
            co_helio = galactic_to_heliocentric(co_dict['rad'],
                                                  co_dict['glat'],
                                                  co_dict['glon'])
            co_dict['px'], co_dict['py'], co_dict['pz'] = co_helio

            ##########
            # Add kicks to NSs and BHs.
            ##########

            # Maxwellian pdf is sqrt(2/pi) * x^2 * exp[-x^2/(2 * a^2)] / a^3,
            # with a scale parameter `a`.
            # The Maxwellian mean is 2 * a * sqrt(2/pi).
            # Here we calculate the scipy.stats `scale` by dividing the
            # user defined mean by the Maxwellian mean.

            NS_idx = np.where(co_dict['rem_id'] == 102)[0]
            NS_kick_speed_scale = NS_kick_speed_mean / (2*np.sqrt(2/np.pi))
            if len(NS_idx) > 0:
                NS_kick_speed = maxwell.rvs(loc=0, scale=NS_kick_speed_scale, size=len(NS_idx))
                NS_kick = utils.sample_spherical(len(NS_idx), NS_kick_speed)
                co_dict['vx'][NS_idx] += NS_kick[0]
                co_dict['vy'][NS_idx] += NS_kick[1]
                co_dict['vz'][NS_idx] += NS_kick[2]

            BH_idx = np.where(co_dict['rem_id'] == 103)[0]
            BH_kick_speed_scale = BH_kick_speed_mean / (2*np.sqrt(2/np.pi))
            if len(BH_idx) > 0:
                BH_kick_speed = maxwell.rvs(loc=0, scale=BH_kick_speed_scale, size=len(BH_idx))
                BH_kick = utils.sample_spherical(len(BH_idx), BH_kick_speed)
                co_dict['vx'][BH_idx] += BH_kick[0]
                co_dict['vy'][BH_idx] += BH_kick[1]
                co_dict['vz'][BH_idx] += BH_kick[2]

            # Add precision to r, b, l
            co_dict['rad'] = utils.add_precision64(co_dict['rad'], -4)
            co_dict['glat'] = utils.add_precision64(co_dict['glat'], -4)
            co_dict['glon'] = utils.add_precision64(co_dict['glon'], -4)

            # Assign vr, mu_b, mu_lcosb.
            co_sph = calc_sph_motion(co_dict['vx'],
                                       co_dict['vy'],
                                       co_dict['vz'],
                                       co_dict['rad'],
                                       co_dict['glat'],
                                       co_dict['glon'])
            co_dict['vr'], co_dict['mu_b'], co_dict['mu_lcosb'] = co_sph

            # Add precision to vr, mu_b, mu_lcosb
            co_dict['vr'] = utils.add_precision64(co_dict['vr'], -4)
            co_dict['mu_b'] = utils.add_precision64(co_dict['mu_b'], -4)
            co_dict['mu_lcosb'] = utils.add_precision64(co_dict['mu_lcosb'], -4)

            # Assign age.
            co_dict['age'] = log_age * np.ones(len(co_dict['vx']))

            ##########
            # Initialize values for compact object photometry.
            # All of the values are np.nan
            # These are all the outputs from the IFMR object in SPISEA.
            ##########
            co_dict['exbv'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_I'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_K'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_J'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_U'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_R'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_B'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_V'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['ubv_H'] = np.full(len(co_dict['vx']), np.nan)

            if additional_photometric_systems is not None:
                if 'ztf' in additional_photometric_systems:
                    co_dict['ztf_g'] = np.full(len(co_dict['vx']), np.nan)
                    co_dict['ztf_r'] = np.full(len(co_dict['vx']), np.nan)
                    co_dict['ztf_i'] = np.full(len(co_dict['vx']), np.nan)

            #########
            # Initialize values for compact object teff, specific gravity and bolometric luminosity
            # All of the values are np.nan
            #########
            co_dict['teff'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['grav'] = np.full(len(co_dict['vx']), np.nan)
            co_dict['mbol'] = np.full(len(co_dict['vx']), np.nan)

            #########
            # Assign feh values to the SPISEA compact objects, drawing randomly from the feh of the Galaxia stars
            #########
            co_dict['feh'] = np.random.choice(star_dict['feh'], size=len(co_dict['vx']), replace=True)

            ##########
            # FIX THE BAD PHOTOMETRY FOR LUMINOUS WHITE DWARFS
            # AND COMPACT OBJECTS WITH COMPANIONS (mags in SPISEA are system magnitudes)
            # For non-dark WDs (the ones from MIST) and compact objects with companions:
            # Approximate extinction from the nearest (in 3-D space) star.
            # Get WD and system photometry from MIST models.
            ##########
            lum_co_sys_idx = np.argwhere(~np.isnan(co_table['m_ubv_I']))

            if len(lum_co_sys_idx) > 0:
                # The extinction to the luminous white dwarfs/luminous systems is calculated
                # by finding the nearest star in the pop_id / age_bin KDTree
                # to the compact object and copying that star's extinction.
                co_xyz = np.array([co_dict['px'][lum_co_sys_idx],
                                     co_dict['py'][lum_co_sys_idx],
                                     co_dict['pz'][lum_co_sys_idx]]).T
                dist, indices = kdt_star_p.query(co_xyz)
                co_dict['exbv'][lum_co_sys_idx] = exbv_arr4kdt[indices.T]

                co_dict['ubv_I'][lum_co_sys_idx] = co_table['m_ubv_I'][lum_co_sys_idx].data
                co_dict['ubv_K'][lum_co_sys_idx] = co_table['m_ukirt_K'][lum_co_sys_idx].data
                co_dict['ubv_J'][lum_co_sys_idx] = co_table['m_ukirt_J'][lum_co_sys_idx].data
                co_dict['ubv_U'][lum_co_sys_idx] = co_table['m_ubv_U'][lum_co_sys_idx].data
                co_dict['ubv_R'][lum_co_sys_idx] = co_table['m_ubv_R'][lum_co_sys_idx].data
                co_dict['ubv_B'][lum_co_sys_idx] = co_table['m_ubv_B'][lum_co_sys_idx].data
                co_dict['ubv_V'][lum_co_sys_idx] = co_table['m_ubv_V'][lum_co_sys_idx].data
                co_dict['ubv_H'][lum_co_sys_idx] = co_table['m_ukirt_H'][lum_co_sys_idx].data
                if additional_photometric_systems is not None:
                    if 'ztf' in additional_photometric_systems:
                        co_dict['ztf_g'][lum_co_sys_idx] = co_table['m_ztf_g'][lum_co_sys_idx].data
                        co_dict['ztf_r'][lum_co_sys_idx] = co_table['m_ztf_r'][lum_co_sys_idx].data
                        co_dict['ztf_i'][lum_co_sys_idx] = co_table['m_ztf_i'][lum_co_sys_idx].data

                # Memory cleaning
                del co_table
                gc.collect()

            # Assign population and object ID.
            co_dict['popid'] = star_dict['popid'][0] * np.ones(len(co_dict['vx']))

            # Multi-Threaded (works in Single-Threaded too)
            with lock:
                co_dict['obj_id'] = np.arange(len(co_dict['vx'])) + next_id_co_val.value
                next_id_co_val.value += len(co_dict['vx'])

    else:
        co_dict = None

    return co_dict, next_id_co_val


def _generate_compound_dtype(obj_arr):
    """
    Create compound datatype by looping over the keys of the obj_arr.
    Assigns integers as datatypes where reasonable, and float64 to the rest
    of the columns

    Parameters
    -----------
    obj_arr : array or None
        Array of stars or compact objects to be binned
        (either co_dict or star_dict)

    Returns
    -------
    compound_dtype : np.dtype
        Numpy datatype with all of the keys of obj_arr
    """

    compound_dtype_arr = []
    for key in obj_arr:
        if key in ['rem_id', 'popid']:  # int16 (up to 32767)
            d = (key, 'i2')
        elif key in ['obj_id']:  # int32 (up to 2147483647)
            d = (key, 'i4')
        else:
            d = (key, 'f8')  # float64
        compound_dtype_arr.append(d)
    compound_dtype = np.dtype(compound_dtype_arr)
    return compound_dtype


def _bin_lb_hdf5(lat_bin_edges, long_bin_edges, obj_arr, output_root, companion_obj_arr = None):
    """
    Given stars and compact objects, sort them into latitude and
    longitude bins. Save each latitude and longitude bin, and the edges that
    define the latitude and longitude bins, as datasets in a compound type hdf5 file.

    Parameters
    -----------
    lat_bin_edges : array
        Edges for the latitude binning (deg)

    long_bin_edges : array
        Edges for the longitude binning (deg)

    obj_arr : array or None
        Array of stars or compact objects to be binned
        (either co_dict or star_dict)

    output_root : str
        The path and name of the hdf5 file,
        without suffix (will be saved as output_root.h5)

    companion_obj_arr : astropy table or None, optional
        Companion table from the ResolvedCluster object. 
        To be used if creating a companion hdf5 file.
        Default None.
        
    Returns
    ------------
    output_root.h5 : hdf5 file
        An compoud type hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects and stars sorted into
        those bins.
    """
    if companion_obj_arr is None:
        # Create compound datatype from obj_arr
        compound_dtype = _generate_compound_dtype(obj_arr)
    else:
        compound_dtype = np.dtype(companion_obj_arr)

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
                                            dtype=compound_dtype)
            else:
                dataset = hf[dset_name]

            ##########
            # Binning the stars and/or compact objects or companions
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
                save_data = np.empty(len(id_lb), dtype=compound_dtype)

                if companion_obj_arr is None:
                    for colname in obj_arr:
                        save_data[colname] = obj_arr[colname][id_lb]
                # If making a companion hd5f file, finds corresponding companions and save them
                else:
                    companion_id_lb = [np.where(companion_obj_arr['system_idx'] == ii)[0] for ii in obj_arr['obj_id'][id_lb]]
                    companion_id_lb = list(np.concatenate(companion_id_lb).ravel()) # Simplifies datastructure
                    if len(companion_id_lb) == 0:
                        continue
                    save_data = np.array(companion_obj_arr[companion_id_lb])

                # Resize the dataset and add data.
                old_size = dataset.shape[0]
                if companion_obj_arr is None:
                    new_size = old_size + len(id_lb)                
                else:
                    new_size = old_size + len(companion_id_lb)                     
                dataset.resize((new_size, ))
                dataset[old_size:new_size] = save_data

            hf.close()

    return

def _no_bins_hdf5(obj_arr, output_root, companion_obj_arr = None):
    """
    Given stars and compact objects save into a compound type hdf5 file.
    No bins, only one dataset (used if doing full sky downsampled)

    Parameters
    -----------
    obj_arr : array or None
        Array of stars or compact objects to be binned
        (either co_dict or star_dict)

    output_root : str
        The path and name of the hdf5 file,
        without suffix (will be saved as output_root.h5)

    companion_obj_arr : astropy table or None, optional
        Companion table from the ResolvedCluster object. 
        To be used if creating a companion hdf5 file.
        Default None.
        
    Returns
    -------
    output_root.h5 : hdf5 file
        An compoud type hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects and stars sorted into
        those bins.
    """
    if companion_obj_arr is None:
        # Create compound datatype from obj_arr
        compound_dtype = _generate_compound_dtype(obj_arr)
    else:
        compound_dtype = np.dtype(companion_obj_arr)

    ##########
    # Loop through the latitude and longitude bins.
    ##########
    # Open our HDF5 file for reading and appending.
    # Create as necessary.
    hf = h5py.File(output_root + '.h5', 'r+')
    # HDF5 dataset name
    dset_name = 'objects'

    # Create data set the first time
    if dset_name not in hf:
        dataset = hf.create_dataset(dset_name, shape=(0,),
                                    chunks=(1e4,),
                                    maxshape=(None,),
                                    dtype=compound_dtype)
    else:
        dataset = hf[dset_name]
    
    ##########
    # Saving the stars and/or compact objects or companions
    ##########
    if obj_arr is not None:
        # Loop over the obj_arr and add all columns into save_data
        save_data = np.empty(len(obj_arr['mass']), dtype=compound_dtype)
        if companion_obj_arr is None:
            for colname in obj_arr:
                save_data[colname] = obj_arr[colname]
        # If making a companion hd5f file, finds corresponding companions and save them
        else:
            save_data = np.array(companion_obj_arr['system_idx'])
        
        # Resize the dataset and add data.
        old_size = dataset.shape[0]
        if companion_obj_arr is None:
            new_size = old_size + len(obj_arr['mass'])                
        else:
            new_size = old_size + len(companion_obj_arr['system_idx'])                     
        dataset.resize((new_size, ))
        dataset[old_size:new_size] = save_data

    hf.close()

    return


def _make_cluster(iso_dir, log_age, currentClusterMass,
                  multiplicity=None, IFMR = 'Raithel18',
                  feh = 0, seed=None, additional_photometric_systems=None):
    """
    Creates SPISEA ResolvedCluster() object.

    Parameters
    -----------
    iso_dir : filepath
        Where are the isochrones stored (for SPISEA)

    log_age : float
        log(age/yr) of the cluster you want to make

    currentClusterMass : float
        Mass of the cluster you want to make (M_sun)

    additional_photometric_systems : list of strs or None, optional
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.

    seed : int or None, optional
         Seed used to sample the kde tree. If set to any number,
         SPISEA will also be forced to use 42 as a
         random seed for calls to ResolvedCluster.
         Default is None.
         
    multiplicity: object or None, optional
        If a resovled multiplicity object is specified, 
        the table will be generated with resolved multiples.
        Default is None.
        
    IFMR : string
        The name of the IFMR object from SPISEA. For more information on these objects see ifmr.py 
        in SPISEA. 
        'Raithel18' = IFMR_Raithel18
        'Spera15' = IFMR_Spera15
        'SukhboldN20' = IFMR_N20_Sukhbold
        
    feh : float
        metallicity [Fe/H] of the cluster you want to make, rounds to nearest value in this list for MIST.
        [-4.0, -3.5, -3.0, -2.5, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.50, -0.25, 0, 0.25, 0.5]
        The IFMR object will run at exactly the specified value.


    Returns
    -------
    cluster : object
        Resolved cluster object from SPISEA.
        
    unmade_cluster_counter_new : int
        Updated number of unmade clusters (<= 100 M_sun)

    unmade_cluster_mass_new: float
        The current mass in the unmade clusters (<= 100 M_sun)

    """
    cluster = None
    
    # Add additional filters to isochrones if additional_photometric_systems
    # contains photometric systems
    my_filt_list = copy.deepcopy(all_filt_list)
    if additional_photometric_systems is not None:
        if 'ztf' in additional_photometric_systems:
            my_filt_list += ['ztf,g', 'ztf,r', 'ztf,i']

    # Calculate the initial cluster mass
    # changed from 0.08 to 0.11 at start because MIST can't handle.
    massLimits = np.array([0.11, 0.5, 120])
    powers = np.array([-1.3, -2.3])
    

    ##########
    # Create the SPISEA table (stars and compact objects).
    #    - it is only sensible to do this for a decent sized cluster.
    ##########
    if currentClusterMass <= 80:
        unmade_cluster_counter = 1
        unmade_cluster_mass = currentClusterMass
                    
    else:
        unmade_cluster_counter = 0
        unmade_cluster_mass = 0
        
        # MAKE isochrone
        # -- arbitrarily chose AKs = 0, distance = 10 pc
        # (irrelevant, photometry not used)
        # Using MIST models to get white dwarfs
        my_iso = synthetic.IsochronePhot(log_age, 0, 10,
                                         evo_model=evolution.MISTv1(),
                                         filters=my_filt_list,
                                         iso_dir=iso_dir,
                                         metallicity=feh)

        # Check that the isochrone has all of the filters in filt_list
        # If not, force recreating the isochrone with recomp=True
        my_iso_filters = [f for f in my_iso.points.colnames if 'm_' in f]
        my_filt_list_fmt = ['m_%s' % f.replace(',', '_') for f in my_filt_list]
        # Checks if the list of filters are different
        if len(set(my_filt_list_fmt) - set(my_iso_filters)) > 0:
            my_iso = synthetic.IsochronePhot(log_age, 0, 10,
                                             evo_model=evolution.MISTv1(),
                                             filters=my_filt_list,
                                             iso_dir=iso_dir,
                                             recomp=True,
                                             metallicity=feh)

        # !!! Keep trunc_kroupa out here !!! Death and destruction otherwise.
        # DON'T MOVE IT OUT!
        trunc_kroupa = imf.IMF_broken_powerlaw(massLimits, powers, multiplicity=multiplicity)

        # MAKE cluster in chunks of currentClusterMass to match 
        # the SPISEA primary star matches to the galaxia mass
        cluster = None
        SPISEA_present_day_star_mass = 0
        while SPISEA_present_day_star_mass < currentClusterMass:
            cluster_chunk_mass = currentClusterMass
            if cluster is None:
                cluster = synthetic.ResolvedCluster(my_iso, trunc_kroupa,
                                            cluster_chunk_mass, ifmr=IFMR_dict[IFMR],
                                            seed=seed, verbose=False)
                _rename_mass_columns_from_spisea(cluster, multiplicity)

            elif cluster is not None:
                cluster_addition = synthetic.ResolvedCluster(my_iso, trunc_kroupa,
                                            cluster_chunk_mass, ifmr=IFMR_dict[IFMR],
                                            seed=seed, verbose=False)
                _rename_mass_columns_from_spisea(cluster_addition, multiplicity)

                if multiplicity is not None:
                    # index of companion assocaited with row of primary, so bump it down the length
                    # of the rest of the cluster
                    cluster_addition.companions['system_idx'] += len(cluster.star_systems)
                    cluster.companions = vstack([cluster.companions, cluster_addition.companions])

                cluster.star_systems = vstack([cluster.star_systems, cluster_addition.star_systems])

            # Note, from here on, mass = current mass, zams_mass = initial mass.
            SPISEA_present_day_star_mass = np.sum(cluster.star_systems[cluster.star_systems['phase'] < 100]['mass'])
            SPISEA_total_mass = np.sum(cluster.star_systems['mass'])
        
        
        # When it overshoots loop through the stars until we get to the appropriate mass
        # since it's in a random order this should be fine
        # this will keep the compact objects that are interspursed along the way, 
        # but their mass is intentionally not counted since we're matching the star mass
        SPISEA_present_day_star_mass_fix = 0
        cluster.star_systems['counting_index'] = np.arange(0, len(cluster.star_systems))
        for idx, object_mass, object_phase in cluster.star_systems.iterrows('counting_index', 'mass', 'phase'):
            if object_phase < 100:
                SPISEA_present_day_star_mass_fix += object_mass
            if SPISEA_present_day_star_mass_fix > currentClusterMass:
                last_index = idx
                break                
        #+1 for difference in index definitions
        cluster.star_systems = cluster.star_systems[:last_index + 1]
        if multiplicity is not None:
            last_multiple_idx = np.where(cluster.star_systems['isMultiple'] == 1)[0][-1]
            # take the last matching index so triples aren't separated
            companion_cutoff_index = np.where(cluster.companions['system_idx'] == last_multiple_idx)[0][-1] + 1
            cluster.companions = cluster.companions[:companion_cutoff_index]
        
        # Only check for clusters greater than 5,000 Msun since
        # 1% of cluster mass could be one star for smaller clusters
        if currentClusterMass > 5*10**4:
            SPISEA_present_day_star_mass_after_matching = np.sum(cluster.star_systems[cluster.star_systems['phase'] < 100]['mass'])
            assert(np.abs((SPISEA_present_day_star_mass_after_matching/currentClusterMass) - 1) < 0.01)

        # Removes companions with unnassigned masses being not massive enough for
        # their isochrone (often high metallicity and old)
        if multiplicity is not None:
            bad_companions = np.where(np.isnan(cluster.companions['mass']) == True)
            bad_companions_primaries = cluster.companions['system_idx'][bad_companions]
            unique_bad_companions_primaries = np.unique(bad_companions_primaries, return_counts = True)
            # Subtract number of companions that are bad
            cluster.star_systems['N_companions'][unique_bad_companions_primaries[0]] -= unique_bad_companions_primaries[1]
            # Set those that have no companions left to non-multiples
            cluster.star_systems['isMultiple'][cluster.star_systems['N_companions'] == 0] = 0
            # Drop companions
            cluster.companions.remove_rows(bad_companions)
        
        
    return cluster, unmade_cluster_counter, unmade_cluster_mass


def _rename_mass_columns_from_spisea(cluster, multiplicity):
    """Rename SPISEA columns:
         SPISEA mass column (initial mass) --> zams_mass
         SPISEA mass_current (current mass) --> mass

    Parameters
    -----------
    cluster : object
        Resolved cluster object from SPISEA.

    multiplicity: object
        If a resovled multiplicity object is specified,
        the table will be generated with resolved multiples.
    """
    cluster.star_systems.rename_column('mass', 'zams_mass')
    cluster.star_systems.rename_column('mass_current', 'mass')
    if multiplicity is not None:
        cluster.companions.rename_column('mass', 'zams_mass')
        cluster.companions.rename_column('mass_current', 'mass')

    return

def _add_multiples(star_zams_masses, cluster, verbose=0):
    """
    Modifies companion table of cluster object to point to Galaxia stars.
    Effectively adds multiple systems with stellar primaries.

    Parameters
    -----------
    star_zams_masses : list
        Galaxia star zams_mass column form the star_dict.

    cluster : object
        Resolved cluster object from SPISEA.

    Ouput:
    -------
    modified_companions : astropy table
        cluster companion table modified to point at Galaxia stars.
        Deleted companions of compact objects and with primary zams masses too large.

    """
    return _add_multiples_some_unmatched(star_zams_masses, cluster, verbose=verbose)

def _add_multiples_some_unmatched(star_zams_masses, cluster, verbose=0):
    """
    Modifies companion table of cluster object to point to Galaxia stars.
    Effectively adds multiple systems with stellar primaries.

    Parameters
    -----------
    star_zams_masses : list
        Galaxia star zams_mass column form the star_dict.

    cluster : object
        Resolved cluster object from SPISEA.

    Returns
    -------
    modified_companions : astropy table
        cluster companion table modified to point at Galaxia stars.
        Deleted companions of compact objects and with primary zams masses too large.

    """
    # FIXME could join table above to improve memory instead of modifying the ss
    # and multiple tables separately

    # Define a timer to keep track of runtimes in this function.
    if verbose > 2:
        t0 = time.time()
        print(f'\t Timer _add_multiples: start {t0:.5f} sec')

    # populate an index column into cluster_ss to preserve original index
    cluster_ss = cluster.star_systems
    cluster_ss['system_idx'] = np.arange(len(cluster_ss))

    # cut down primaries table with the following conditions:
    #  - only multiples
    #  - only stellar primaries (no dark compact objects)
    #  - only spisea stars less massive than the most massive Galaxia star.
    # These conditions are applied successively.
    cond_multiple = cluster_ss['isMultiple'] == True
    cond_stellar_prim = cond_multiple & (cluster_ss['phase'] < 100)
    cond_not_too_massive = cond_stellar_prim & (cluster_ss['zams_mass'] <= np.max(star_zams_masses))
    if verbose > 3:
        print(f'_add_multiples: Trimming from total of {len(cluster_ss)} SPISEA primaries.')
        print(f'  {np.sum(cond_multiple)} are multiple')
        print(f'  {np.sum(cond_stellar_prim)} are multiple and stars')
        print(f'  {np.sum(cond_not_too_massive)} are multiple, stars, and not too massive')
    cluster_ss = cluster_ss[cond_not_too_massive]
    
    
    # For each primary, find the closest, higher-mass Galaxia star.
    if verbose > 2:
        print(f'\t Timer _add_multiples: test1 {time.time() - t0:.5f} sec')
    closest_index_arr, closest_mass_diff = match_companions(star_zams_masses,
                                                            cluster_ss['zams_mass'].data,
                                                            verbose=verbose)
    if verbose > 2:
        print(f'\t Timer _add_multiples: test2 {time.time() - t0:.5f} sec')

    # Add columns for the matched galaxia ID and primary zams mass difference.
    cluster_ss['galaxia_id'] = closest_index_arr
    cluster_ss['zams_mass_match_diff'] = closest_mass_diff

    # Trim out the non-matches
    cond_matched = closest_index_arr != -1
    if verbose > 3:
        print(f'_add_multiples: Trimming from total of {len(cluster_ss)} SPISEA primaries in matching.')
        print(f'  {np.sum(cond_matched)} are matched')
    cluster_ss = cluster_ss[cond_matched]

    # join tables to keep only those companions that pass our cuts.
    companion_tmp_df = cluster.companions.to_pandas().set_index(['system_idx'])
    event_table_df = cluster_ss.to_pandas().set_index(['system_idx'])
    companion_tmp_df_joined = companion_tmp_df.join(event_table_df, lsuffix='', rsuffix='_prim', how='inner')
    
    # for cc in range(len(companion_tmp_df_joined)):
    #     print(cc)
    #     print(companion_tmp_df_joined['zams_mass_match_diff'].values[cc])
    #     print(star_zams_masses[companion_tmp_df_joined['galaxia_id'].values[cc]])
    #     print(companion_tmp_df_joined['zams_mass_prim'].values[cc])
    #     print(star_zams_masses[companion_tmp_df_joined['galaxia_id'].values[cc]] - companion_tmp_df_joined['zams_mass_prim'].values[cc])
    #     np.testing.assert_almost_equal(companion_tmp_df_joined['zams_mass_match_diff'].values[cc],
    #                                star_zams_masses[companion_tmp_df_joined['galaxia_id'].values[cc]] - companion_tmp_df_joined['zams_mass_prim'].values[cc],
    #                                decimal=2)

    # defines the companion Table as the new cut down table
    # reset index adds system_idx back as its own column
    # doing list(cluster.companions.columns) means that it only takes columns that were in cluster.companions
    # all values that were nan were switched to None in pandas, so .filled(np.nan) switches them back to nan.
    final_columns = list(cluster.companions.columns) + ['galaxia_id', 'zams_mass_match_diff', 'zams_mass_prim']
    modified_companions = Table.from_pandas(
        companion_tmp_df_joined.reset_index()[final_columns]).filled(np.nan)
    
    # confirm all compnions were assigned to a popsycle primary
    assert np.sum(np.isnan(modified_companions['galaxia_id'])) == 0
    assert np.sum(modified_companions['galaxia_id'] < 0) == 0

    # Reorganize the columns
    modified_companions['spisea_idx'] = modified_companions['system_idx']
    #modified_companions['spisea_idx'] = copy.deepcopy(modified_companions['system_idx'])
    modified_companions['system_idx'] = modified_companions['galaxia_id'].astype(int)
    del modified_companions['galaxia_id']
    #del cluster.star_systems['system_idx']

    if modified_companions is None:
        return None

    return modified_companions

def _add_multiples_all(star_zams_masses, cluster, verbose=0):
    """
    Modifies companion table of cluster object to point to Galaxia stars.
    Effectively adds multiple systems with stellar primaries.

    Parameters
    -----------
    star_zams_masses : list
        Galaxia star zams_mass column form the star_dict.

    cluster : object
        Resolved cluster object from SPISEA.

    Returns
    -------
    modified_companions : astropy table
        cluster companion table modified to point at Galaxia stars.
        Deleted companions of compact objects and with primary zams masses too large.

    """
    # FIXME could join table above to improve memory instead of modifying the ss
    # and multiple tables separately

    # Define a timer to keep track of runtimes in this function.
    t0 = time.time()
    if verbose > 2:
        print(f'\t Timer _add_multiples: start {t0:.5f} sec')

    # populate an index column into cluster_ss to preserve original index
    cluster_ss = cluster.star_systems
    cluster_ss['system_idx'] = np.arange(len(cluster_ss))

    # cut down table to only multiples
    cond_multiple = cluster_ss['isMultiple'] == True
    cluster_ss = cluster_ss[cond_multiple]
    
    # join tables to cut on companions table
    companion_tmp_df = cluster.companions.to_pandas().set_index(['system_idx'])
    event_table_df = cluster_ss.to_pandas().set_index(['system_idx'])
    companion_tmp_df_joined = companion_tmp_df.join(event_table_df, lsuffix='', rsuffix='_prim', how='inner')

    # cut down table to only stellar primaries
    companion_tmp_df_joined = companion_tmp_df_joined[companion_tmp_df_joined['phase_prim'] < 100]
    cluster_ss = cluster_ss[cluster_ss['phase'] < 100]

    # ALL SPISEA SYSTEMS WITH PRIMARIES MORE MASSIVE THAN THE MOST MASSIVE GALAXIA
    # SYSTEM ARE DROPPED!!!!!!!!! (0-2% of systems).
    companion_tmp_df_joined = companion_tmp_df_joined[companion_tmp_df_joined['zams_mass_prim'] <= np.max(star_zams_masses)]
    cluster_ss = cluster_ss[cluster_ss['zams_mass'] <= np.max(star_zams_masses)]
    too_big = np.sum(cluster_ss['zams_mass'] > np.max(star_zams_masses))
    
    # defines the companion Table as the new cut down table
    # reset index adds system_idx back as its own column
    # doing list(cluster.companions.columns) means that it only takes columns that were in cluster.companions
    # all values that were nan were switched to None in pandas, so .filled(np.nan) switches them back to nan.
    modified_companions = Table.from_pandas(companion_tmp_df_joined.reset_index()[list(cluster.companions.columns)]).filled(np.nan)

    # Place new system_idx into temporary column initiated with NaN
    modified_companions['system_idx_tmp'] = np.nan
    modified_companions['zams_mass_match_diff'] = np.nan

    if verbose > 2:
        print(f'\t Timer _add_multiples: test1 {time.time() - t0:.5f} sec')
    closest_index_arr, mass_diff = _match_companions_kdtree(star_zams_masses, cluster_ss['zams_mass'])
    if verbose > 2:
        print(f'\t Timer _add_multiples: test2 {time.time() - t0:.5f} sec')

    # loop through each of the SPISEA cluster and match them to a galaxia star
    for ii, index in enumerate(cluster_ss['system_idx']):
        closest_index = closest_index_arr[ii]
        companion_indicies = np.where(modified_companions['system_idx'] == index)[0]
        # points companions to nearest-in-mass Galaxia primary to SPISEA primary
        for jj in companion_indicies:
            modified_companions[jj]['system_idx_tmp'] = closest_index
            modified_companions[jj]['zams_mass_match_diff'] = star_zams_masses[closest_index] - cluster_ss[ii]['zams_mass']
    
    # confirm all compnions were assigned to a popsycle primary
    assert np.sum(np.isnan(modified_companions['system_idx_tmp'])) == 0
    assert len(set(modified_companions['system_idx_tmp'])) == len(closest_index_arr)
    modified_companions['system_idx'] = modified_companions['system_idx_tmp'].astype(int)
    del modified_companions['system_idx_tmp']
    del cluster.star_systems['system_idx']

    if modified_companions is None:
        return None

    num_companions = len(cluster.companions)
    if verbose > 1:
        print("Total companions, too big, too big fraction:", num_companions,
              too_big, too_big / num_companions)

    return modified_companions

def match_companions(star_zams_masses, SPISEA_primary_zams_masses, verbose=0):
    """
    Matches galaxia stars and SPISEA stellar primaries by closest match.
    Note this is after all SPISEA stars with mass greater than the largest
    galaxia mass have been dropped.

    Parameters
    ----------
    star_zams_masses : list
        Galaxia star zams mass column form the star_dict.

    SPISEA_primary_zams_masses : object
        Astropy column 'zams_mass' from the cluster.star_systems cut down to stellar
        primaries with companions less massive than the most massive galaxia star.

    Returns
    -------
    closest_index_arr : array
        Cloest index in the galaxia mz_galaxia (used to index into mz_galaxia)
        with the length of mz_spisea. Note that an index of -1 indicates this star
        cannot be matched within 20% of the desired mass.
    mass_diff : array
        The difference between the galaxia star mass and the spisea primary mass
        (m_galaxia - m_spisea) for every input spisea star.  Note that a value of np.nan
        is used for unmatched stars.
    """

    closest_index_arr, mass_diff = _match_companions_kdtree_nonneg(star_zams_masses,
                                                                   SPISEA_primary_zams_masses,
                                                                   verbose=verbose)

    return closest_index_arr, mass_diff

def _match_companions_kdtree(star_zams_masses, SPISEA_primary_zams_masses, verbose=0):
    """
    Matches galaxia stars and SPISEA stellar primaries by closest match.
    Note this is after all SPISEA stars with mass greater than the largest
    galaxia mass have been dropped.

    This algorithm utilizes a KDTree with a euclidean distance on mass.
    It is the fastest approach and all stars will return masses. However, it
    will return matches such that m_galaxia < m_spisea, which is not correct.
    Matches may also have very large mass differences.  Duplicates are not
    allowed and the search continues until all matches are unduplicated.

    Parameters
    -----------
    star_zams_masses : list
        Galaxia star zams mass column form the star_dict.

    SPISEA_primary_zams_masses : object
        Astropy column 'zams_mass' from the cluster.star_systems cut down to stellar
        primaries with companions less massive than the most massive galaxia star.

    Returns
    -------
    closest_index_arr : array
        Cloest index in the galaxia mz_galaxia (used to index into mz_galaxia)
        with the length of mz_spisea. Note that an index of -1 indicates this star
        cannot be matched within 20% of the desired mass.
    mass_diff : array
        The difference between the galaxia star mass and the spisea primary mass
        (m_galaxia - m_spisea) for every input spisea star.  Note that a value of np.nan
        is used for unmatched stars.  This function does not return mass_diff (=None).
    """

    # prepare KDTree of Galaxia masses for mass matching
    star_mass_tree = np.expand_dims(star_zams_masses, axis=1)
    galaxia_mass_tree = cKDTree(star_mass_tree)
    
    # search the tree with SPISEA masses for their nearest match
    cluster_search_mass = np.expand_dims(SPISEA_primary_zams_masses, axis=1)
    k = 1
    _, closest_index_arr = galaxia_mass_tree.query(cluster_search_mass, k=k)
    print('Starting with %i SPISEA to Galaxia mass matches' % len(closest_index_arr))

    # find the number of matches that are duplicates.
    # indexes returns the 1st star_system index that
    # has the "unique" closest_index. This is not
    # necessarily the closest match.
    _, indexes, counts = np.unique(closest_index_arr,
                                   return_index=True,
                                   return_counts=True)
    cond = counts > 1
    print('-- Found %i duplicates at k = %i' % (np.sum(cond), k))

    # grab indicies where the first duplicate is located
    nonunique_indicies = indexes[cond]

    # while there are duplicates...
    while np.sum(cond) > 0:
        # increase the search to one neighbor further away and
        # only search on those masses that were duplicates
        k += 1
        _, next_closest_index_arr = galaxia_mass_tree.query(cluster_search_mass[nonunique_indicies], k=k)

        # extract the neighbor that is furthest away within k-neighbors
        next_closest_index_arr = [i[-1] for i in next_closest_index_arr]

        # assign that value to the duplicates in the closest_index_arr
        closest_index_arr[nonunique_indicies] = next_closest_index_arr

        # count the number of duplicates that remain
        _, indexes, counts = np.unique(closest_index_arr,
                                       return_index=True,
                                       return_counts=True)
        cond = counts > 1
        #print('-- Found %i duplicates at k = %i' % (np.sum(cond), k))
        nonunique_indicies = indexes[cond]
    
        
    return closest_index_arr, None

def _match_companions_diff_array(star_zams_masses, SPISEA_primary_zams_masses, verbose=0):
    """
    Matches galaxia stars and SPISEA stellar primaries by closest match.
    Note this is after all SPISEA stars with mass greater than the largest
    galaxia mass have been dropped.

    This algorithm utilizes a 2D array to calculate mass differences for every pair
    of stars. It also only allows matching for stars where m_galaxia >= m_spisea.
    Duplicates are not allowed and the search continues until all matches are unduplicated.
    The algorithm is somewhat slower and consumes too much memory for our purposes.

    Parameters
    -----------
    star_zams_masses : list
        Galaxia star zams mass column form the star_dict.

    SPISEA_primary_zams_masses : object
        Astropy column 'zams_mass' from the cluster.star_systems cut down to stellar
        primaries with companions less massive than the most massive galaxia star.

    Returns
    -------
    closest_index_arr : array
        Cloest index in the galaxia mz_galaxia (used to index into mz_galaxia)
        with the length of mz_spisea. Note that an index of -1 indicates this star
        cannot be matched within 20% of the desired mass.

    mass_diff : array
        The difference between the galaxia star mass and the spisea primary mass
        (m_galaxia - m_spisea) for every input spisea star.  Note that a value of np.nan
        is used for unmatched stars.
    """
    if verbose > 3:
        print('memory test 0')
    
    # Shorter variable names
    m_g = np.asarray(star_zams_masses)
    m_s = np.asarray(SPISEA_primary_zams_masses)
    
    if verbose > 3:
        print('memory test 1')
    
    # Pre-sort both.
    sdx_s = np.argsort(m_s)  # connection between sorted and orig position
    sdx_g = np.argsort(m_g)
    m_s_sorted = m_s[sdx_s]
    m_g_sorted = m_g[sdx_g]
    
    if verbose > 3:
        print('memory test 2')

    # Make a big 2D array of all the differences.
    # Shape is [N_g, N_s]
    diff = m_g_sorted[:, np.newaxis] - m_s_sorted
    
    if verbose > 3:
        print('memory test 3')

    # Append a 1000 Msun diff row so we ALWAYS have a match.
    # Matching to the last row means "no match".
    dummy_diffs = np.array([[1000] * len(m_s), ])
    diff = np.append(diff, dummy_diffs, axis=0)
    idx_nomatch = diff.shape[0] - 1  # index of the dummy case
    
    if verbose > 3:
        print('memory test 4')

    # Only allow matches where the Galaxia mass is larger than the SPISEA mass.
    diff[diff < 0] = np.nan # Set zeros (lower triangle) to nans.
    closest_index_arr_sorted = np.nanargmin(diff, axis=0)
    assert closest_index_arr_sorted.shape[0] == len(m_s)
    del diff
    
    if verbose > 3:
        print('memory test 5')

    # Now loop through and make sure we only have unique matches.
    # Essentially, if two spisea stars match to the same Galaxia object,
    # we will take the higher-mass spisea star and walk up in mass in the Galaxia
    # table until we find a unique match.
    k = 0
    uni, indexes, counts = np.unique(closest_index_arr_sorted,
                                   return_index=True,
                                   return_counts=True)
    
    if verbose > 3:
        print('memory test 6')

    # Count up the number of duplicates, not counting "no-matches"
    cond = (counts > 1) & (uni != idx_nomatch)
    n_dups = np.sum(cond)
    n_nomatch = np.sum(uni == idx_nomatch)
    if verbose > 3:
        print(f'match_companions: Found {n_dups} duplicates in round {k}')
        print(f'match_companions: Found {n_nomatch} with no matches in round {k}')

    nonunique_indices = indexes[cond]

    # loop through until all duplicates are matched uniquely or we run out
    # of galaxia stars at the highest mass. These become "no matches".
    while n_dups > 0:
        # increment up to the next most-massive galaxia star for the duplicates.
        k += 1
        closest_index_arr_sorted[nonunique_indices] += 1

        # Update duplicates check.
        uni, indexes, counts = np.unique(closest_index_arr_sorted,
                                       return_index=True,
                                       return_counts=True)

        cond = (counts > 1) & (uni != idx_nomatch)
        n_dups = np.sum(cond)
        n_nomatch = np.sum(uni == idx_nomatch)
        if verbose > 3:
            print(f'match_companions: Found {n_dups} duplicates in round {k}')
            print(f'match_companions: Found {n_nomatch} with no matches in round {k}')

        nonunique_indices = indexes[cond]


    # Fix up the non-matches.
    good_matches = np.argwhere(closest_index_arr_sorted != idx_nomatch)

    # The final indices are still on the sorted array. Reverse to index on original input arrays.
    # Note that index = -1 means no match.
    closest_index_arr = np.ones_like(sdx_s) * -1
    closest_index_arr[sdx_s[good_matches]] = sdx_g[closest_index_arr_sorted[good_matches]]

    mass_diff = np.full_like(m_s, np.nan)
    mass_diff[sdx_s[good_matches]] = m_g_sorted[closest_index_arr_sorted[good_matches]] - m_s_sorted[good_matches]

    return closest_index_arr, mass_diff

def _match_companions_kdtree_nonneg(mz_galaxia_in, mz_spisea_in, max_frac_mass_diff=0.2, verbose=0):
    """
    Matches galaxia stars and SPISEA stellar primaries by closest match.
    Note this is after all SPISEA stars with mass greater than the largest
    galaxia mass have been dropped.

    This algorithm utilizes a KDTree for matching; but then properly disallows
    the case where m_galaxia < m_spisea. Duplicates are not allowed and the search
    continues until all matches are unduplicated. However, there is a stopping
    criteria set for the maximum allowed fractional mass difference:
        (m_galaxia - m_spisea) / m_galaxia
    Once all remaining unmatched stars have a fractiona mass difference value
    that exceeds 0.2 (default), the search is stopped and all remaining stars are
    declared unmatched. This can amount to ~10% of stars.
    
    Parameters
    -----------
    mz_galaxia_in : list
        Galaxia star zams mass column form the star_dict.

    mz_spisea_in : object
        Astropy column 'zams_mass' from the cluster.star_systems cut down to stellar
        primaries with companions less massive than the most massive galaxia star.

    Returns
    -------
    closest_index_arr : array
        Cloest index in the galaxia mz_galaxia (used to index into mz_galaxia)
        with the length of mz_spisea. Note that an index of -1 indicates this star
        cannot be matched within 20% of the desired mass.

    mass_diff : array
        The difference between the galaxia star mass and the spisea primary mass
        (m_galaxia - m_spisea) for every input spisea star.  Note that a value of np.nan
        is used for unmatched stars.
    
    """
    # prepare KDTree of Galaxia masses for mass matching. KDtree requires 2d arrays (N_stars, 1)
    mz_galaxia = np.expand_dims(mz_galaxia_in, axis=1)
    mz_spisea = np.expand_dims(mz_spisea_in, axis=1)
    galaxia_mass_tree = cKDTree(mz_galaxia)
    if verbose > 0:
        print(f'match_companions: {mz_galaxia.shape[0]} galaxia masses, {mz_spisea.shape[0]} spisea masses')

    # Setup our output variables.
    # Note: closest_index_arr can take the following values:
    #    -2     not yet declared matched or unmatched
    #    -1     unmatched and not to be searched further (e.g. dm/m > 0.2)
    #    >= 0   matched
    closest_index_arr = np.full(mz_spisea.shape[0], -2, dtype=int)
    mass_diff = np.full(mz_spisea.shape[0], np.nan, dtype=np.float16)

    # Start with everything unmatched
    unmatched_indices = np.arange(mz_spisea.shape[0])
    k = 1  # Start with the first nearest neighbor.

    # while there are un-matched stars, keep going.
    while len(unmatched_indices) > 0 and k < mz_galaxia.shape[0]:
        if verbose > 3: print(f'match_companions: k = {k}')

        # Search the KDtree for matches
        _, new_closest_index_arr = galaxia_mass_tree.query(mz_spisea[unmatched_indices], k=k)
        
        # extract the neighbor that is furthest away within k-neighbors
        if len(new_closest_index_arr.shape) > 1:
            new_closest_index_arr = np.array([i[-1] for i in new_closest_index_arr])

        # *** Condition *** Enforce positive mass difference.
        # Calculate the current mass difference.
        new_mass_diff = mz_galaxia[new_closest_index_arr, 0] - mz_spisea[unmatched_indices, 0]

        # Identify the good matches that have m_galaxia < m_spisea.
        new_matches = np.where(new_mass_diff >= 0)[0]
        if verbose > 3: print(f'-- Found {len(unmatched_indices) - len(new_matches)} m_galaxia < m_spisea')

        # For new, good matches, assign the closest index value and mass difference.
        # Otherwise, index = -1 and mass_diff = nan to indicate negative mass diff stars.
        closest_index_arr[unmatched_indices[new_matches]] = new_closest_index_arr[new_matches]
        mass_diff[unmatched_indices[new_matches]] = new_mass_diff[new_matches]

        # *** Condition *** Enforce mass-matching within dm/m < 0.2.
        new_frac_mass_diff = np.abs(new_mass_diff / mz_galaxia[new_closest_index_arr, 0])
        bad_matches = np.where(new_frac_mass_diff > max_frac_mass_diff)[0]
        closest_index_arr[unmatched_indices[bad_matches]] = -1   # Declare these unmatched forever.
        mass_diff[unmatched_indices[bad_matches]] = np.nan
        if verbose > 3:
            print(f'-- Found {len(bad_matches)} new matches with dm/m > {max_frac_mass_diff} to call unmatched.')
            print(f'-- Fractional abs mass diff range = {new_frac_mass_diff.min():.3f} - {new_frac_mass_diff.max():.2f}')

        # *** Condition *** Enforce no duplicates.
        # count the number of duplicates that remain in the whole array.
        uni_vals, indexes, idx_inverse, counts = np.unique(closest_index_arr,
                                              return_index=True,
                                              return_inverse=True,
                                              return_counts=True)
        
        cond_dup = np.logical_and(counts > 1, uni_vals >= 0)  # Only consider dupes among good matches.
        if verbose > 3: print(f'-- Found {np.sum(cond_dup)} duplicates')
    
        # grab indices where the first duplicate is located
        nonunique_indices = indexes[cond_dup]

        # reset the duplicates to unmatched (need for stopping condition below).
        closest_index_arr[nonunique_indices] = -2
        mass_diff[nonunique_indices] = np.nan

        # Fix cases where there were duplicates with >2 counts.
        tdx = np.where((counts > 2) & (uni_vals >= 0))[0]
        if verbose > 3: print(f'-- Found {len(tdx)} duplicates with count > 2')
        for tt in tdx:
            cdx = np.where(idx_inverse == tt)[0]

            # print(f'tt={tt}: cnt={counts[tt]} idx={indexes[tt]}')
            # print(f'         cdx={cdx}')

            # We have already corrected one of them. So for count=3, there should be 2 left.
            # For count=4, there should be 3 left. But we want to keep one of the matches around.
            # Arbitrarily choose the last one to keep around.
            if len(cdx) > 1:
                closest_index_arr[cdx[:-1]] = -2
                mass_diff[cdx[:-1]] = np.nan

        # # temporary to verify that no duplicates remain (except for unmatched -1).
        # uni_vals2, indexes2, idx_inverse2, counts2 = np.unique(closest_index_arr,
        #                                                    return_index=True,
        #                                                    return_inverse=True,
        #                                                    return_counts=True)
        # pdb.set_trace()
        # assert np.max(counts2[1:]) == 1

        # grab indices of all stars that remain to be matched
        # indicated by -2 closest index.
        unmatched_indices = np.where(closest_index_arr == -2)[0]
        if verbose > 3: print(f'-- Found {len(unmatched_indices)} still unmatched after this round')

        # increase the search to one neighbor further away for next round.
        k += 1

    if verbose > 0:
        N_total = len(closest_index_arr)
        N_matched = np.sum(closest_index_arr >= 0)
        #N_unmatched = np.sum(closest_index_arr == -1)
        print(f'match_companions: Maximum k = {k-1} with {N_matched} of {N_total} matched.')

    return closest_index_arr, mass_diff


def _make_companions_table(cluster, star_dict, co_dict,
                           additional_photometric_systems = None, t0 = 0,
                           verbose=0):
    """
    Makes companions table by:
        1. Taking companions of compact objects and pointing their 
        system_idx to the co_dict obj_id (so not dependent on position in list)
        2. Matching SPISEA stellar primaries to galaxia primaries by mass
        and pointing the companion system_idx to the galaxia obj_id
        3. Combining those two tables into companions_table.

    Also modifies the following star_dict columns:
        1. systemMass = primary mass + companion masses (current mass)
        2. Sets isMultiple = 1 for all systems with companions
        3. Sets N_companions column based on duplication of system_idx
        in companions table
        4. Changes magnitudes to system magnitudes

    Parameters
    -----------
    cluster : object
        Resolved cluster object from SPISEA.
    
    star_dict : dictionary
        The number of entries for each key is the number of (primary) stars.
    
    co_dict : dictionary
        Keys are the same as star_dict, just for compact objects.

    additional_photometric_systems : list of strs or None, optional
        The name of the photometric systems which should be calculated from
        Galaxia / SPISEA's ubv photometry and appended to the output files.
        
    t0 : float, optional
        Initial time for timing purposes. Default is 0.

    verbose : int, optional
        Print out more verbose statements for higher numbers. Coarse timing reported
        with verbose=2. Fine timing reported for verbose=4. Default is 0.

    Returns
    -------
    star_dict : dictionary
        The number of entries for each key is the number of (primary) stars.
        Modified systemMass, isMultiple, and N_companion columns.
        
    companions_table : astropy table
        Companions table with columns from SPISEA, but points system_idx to
        galaxia obj_id.

    """
    if cluster:
        if verbose > 3: print(f'test1 {time.time() - t0:.2f} sec')

        ###########
        # Treats compact object companions
        ###########
        
        # Makes a separate companion table with primaries as compact objects
        # Points the system_idx to obj_id instead of idx
        if co_dict is not None and cluster.companions:
            
            # Build a list of systems for every companion (will be repeat systems). Shape = N_companions
            # We will call this the "dup_sys" table for duplicate systems (duplicated x N_companions)
            clust_dup_sys = cluster.star_systems[cluster.companions['system_idx']]
            # Filter the dup_sys table to just those with CO primaries.
            co_idx_in_dup_sys_table = np.where(clust_dup_sys['phase'] > 100)[0]
            compact_companions = cluster.companions[co_idx_in_dup_sys_table]

            # Make a new co_dict table that still preserves N_companions.
            co_idx_in_sys_table = np.where(cluster.star_systems['phase'] > 100)[0]
            co_dict_tmp = cluster.star_systems[co_idx_in_sys_table]

            # Repeat the indices into co_dict times the number of companions
            co_idx = np.arange(0, len(co_dict_tmp))
            co_idx_for_dup_sys_table = np.repeat(co_idx, co_dict_tmp['N_companions'])

            # Reset the companion system_idx to be the correct obj_id
            compact_companions['system_idx'] = co_dict['obj_id'][co_idx_for_dup_sys_table]

            # sums mass of companions if there are triples
            grouped_co_companions = compact_companions.group_by(['system_idx'])
            CO_companions_system_mass = grouped_co_companions['mass'].groups.aggregate(np.sum)
            grouped_system_idxs = np.array(grouped_co_companions.groups.keys['system_idx'])
            # Returns the intersecting obj_ids/system_idxs, indices in CO_table that correspond with the overlap, 
            # and indices in companion_table that overlap (this last one should be just np.arange(len(companion_table)))
            CO_idx_w_companions = np.intersect1d(np.array(co_dict['obj_id']), grouped_system_idxs, return_indices = True)
            
            # indexes on grouped_companions.groups.keys since this is one system_idx per system 
            # (i.e. no duplicated if there are triples)
            co_dict['systemMass'][CO_idx_w_companions[1]] += CO_companions_system_mass

            if verbose > 3: print(f'test2 {time.time() - t0:.2f} sec')

            del co_dict_tmp
        
        ###########
        # Treats stellar companions and joins stellar and CO companions
        ###########
        # First matches galaxia and SPISEA primaries
        companions_table = _add_multiples(star_zams_masses=star_dict['zams_mass'], cluster=cluster)
        if verbose > 3: print(f'test3 {time.time() - t0:.2f} sec')
        if companions_table:
            # sums mass of companions if there are triples
            grouped_companions = companions_table.group_by(['system_idx'])
            companions_system_mass = grouped_companions['mass'].groups.aggregate(np.sum)
            group_companions_system_idxs = grouped_companions.groups.keys['system_idx']
            # indexes on grouped_companions.groups.keys since this is one system_idx per system 
            # (i.e. no duplicated if there are triples)
            star_dict['systemMass'][group_companions_system_idxs] += companions_system_mass
            star_dict['isMultiple'][group_companions_system_idxs] = 1

            # Adds N companions column to star_dict
            # this column already existed in the main cluster
            # but remaking due to indexing changes
            N_companions_dictionary = Counter(companions_table['system_idx'])
            star_dict['N_companions'][companions_table['system_idx']] = itemgetter(*companions_table['system_idx'])(N_companions_dictionary)
            
            # Changes luminosities to be system luminosities
            companions_system_m_ubv_I = grouped_companions['m_ubv_I'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_K = grouped_companions['m_ukirt_K'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_J = grouped_companions['m_ukirt_J'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_U = grouped_companions['m_ubv_U'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_R = grouped_companions['m_ubv_R'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_B = grouped_companions['m_ubv_B'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_V = grouped_companions['m_ubv_V'].groups.aggregate(binary_utils.add_magnitudes)
            companions_system_m_ubv_H = grouped_companions['m_ukirt_H'].groups.aggregate(binary_utils.add_magnitudes)
            if additional_photometric_systems is not None:
                if 'ztf' in additional_photometric_systems:
                    companions_system_m_ztf_g = grouped_companions['m_ztf_g'].groups.aggregate(binary_utils.add_magnitudes)
                    companions_system_m_ztf_r = grouped_companions['m_ztf_r'].groups.aggregate(binary_utils.add_magnitudes)
                    companions_system_m_ztf_i = grouped_companions['m_ztf_i'].groups.aggregate(binary_utils.add_magnitudes)
                    
            
            star_dict['ubv_I'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_I'][group_companions_system_idxs], companions_system_m_ubv_I])
            star_dict['ubv_K'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_K'][group_companions_system_idxs], companions_system_m_ubv_K])
            star_dict['ubv_J'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_J'][group_companions_system_idxs], companions_system_m_ubv_J])
            star_dict['ubv_U'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_U'][group_companions_system_idxs], companions_system_m_ubv_U])
            star_dict['ubv_R'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_R'][group_companions_system_idxs], companions_system_m_ubv_R])
            star_dict['ubv_B'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_B'][group_companions_system_idxs], companions_system_m_ubv_B])
            star_dict['ubv_V'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_V'][group_companions_system_idxs], companions_system_m_ubv_V])
            star_dict['ubv_H'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ubv_H'][group_companions_system_idxs], companions_system_m_ubv_H])
            if additional_photometric_systems is not None:
                if 'ztf' in additional_photometric_systems:
                    star_dict['ztf_g'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ztf_g'][group_companions_system_idxs], companions_system_m_ztf_g])
                    star_dict['ztf_r'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ztf_r'][group_companions_system_idxs], companions_system_m_ztf_r])
                    star_dict['ztf_i'][group_companions_system_idxs] = binary_utils.add_magnitudes([star_dict['ztf_i'][group_companions_system_idxs], companions_system_m_ztf_i])
            
            # Switch companion table to point to obj_id instead of idx
            companions_table['system_idx'] = star_dict['obj_id'][companions_table['system_idx']]

            # Makes companions with stellar and compact primaries into one table
            if co_dict is not None:
                companions_table = (vstack([companions_table, compact_companions], join_type='inner'))

            if verbose > 3: print(f'test4 {time.time() - t0:.2f} sec')
    else:
        companions_table = None
            
    return star_dict, companions_table
            
        

    
    

############################################################################
########### Candidate event calculation and associated functions ###########
############################################################################

def _check_calc_events(hdf5_file, output_root2,
                       radius_cut, obs_time, n_obs, theta_frac,
                       blend_rad, n_proc, overwrite, hdf5_file_comp):
    """
    Checks that the inputs of calc_events are valid

    Parameters
    -----------
    hdf5_file : str
        Name of the HDF5 file.

    output_root2 : str
        The name for the h5 file

    radius_cut : float
        Initial radius cut, in ARCSECONDS.

    obs_time : float
        Survey duration, in DAYS.

    n_obs : int
        Number of observations.

    theta_frac : float
        Another cut, in multiples of Einstein radii.

    blend_rad : float
        Stars within this distance of the lens are said to be blended.
        Units are in ARCSECONDS.

    n_proc : int
        Number of processors to use. Should not exceed the number of cores.
        Default is one processor (no parallelization).

    overwrite : bool
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.
        Default is False.
        
    hdf5_file_comp: str
        String of hdf5 file of companion events created in perform_pop_syn().
        Default is None.
    """

    if not isinstance(hdf5_file, str):
        raise Exception('hdf5_file (%s) must be a string.' % str(hdf5_file))

    if hdf5_file[-3:] != '.h5':
        raise Exception('hdf5_file (%s) must be an h5 file.' % str(hdf5_file))

    if not isinstance(output_root2, str):
        raise Exception('output_root2 (%s) must be a string.' % str(output_root2))

    if not isinstance(radius_cut, int):
        if not isinstance(radius_cut, float):
            raise Exception('radius_cut (%s) must be an integer or a float.' % str(radius_cut))

    if not isinstance(obs_time, int):
        if not isinstance(obs_time, float):
            raise Exception('obs_time (%s) must be an integer or a float.' % str(obs_time))

    if not isinstance(blend_rad, int):
        if not isinstance(blend_rad, float):
            raise Exception('blend_rad (%s) must be an integer or a float.' % str(blend_rad))

    if not isinstance(n_obs, int):
        raise Exception('n_obs (%s) must be an integer.' % str(n_obs))

    if not isinstance(n_proc, int):
        raise Exception('n_proc (%s) must be an integer.' % str(n_proc))

    if not isinstance(overwrite, bool):
        raise Exception('overwrite (%s) must be a boolean.' % str(overwrite))

    if not isinstance(theta_frac, int):
        if not isinstance(theta_frac, float):
            raise Exception('theta_frac (%s) must be an integer or a float.' % str(theta_frac))

    if not isinstance(hdf5_file_comp, str):
        if not isinstance(hdf5_file_comp, type(None)):
            raise Exception('hdf5_file_comp (%s) must be a str or a NoneType.' % str(hdf5_file_comp))


def calc_events(hdf5_file, output_root2,
                radius_cut=2, obs_time=1000, n_obs=101, theta_frac=2,
                blend_rad=0.65, n_proc=1,
                overwrite=False, hdf5_file_comp=None):
    """
    Calculate microlensing events

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file.

    output_root2 : str
        The name for the h5 file

    radius_cut : float
        Initial radius cut, in ARCSECONDS.

    obs_time : float
        Survey duration, in DAYS.

    n_obs : int
        Number of observations.

    theta_frac : float
        Another cut, in multiples of Einstein radii.

    blend_rad : float
        Stars within this distance of the lens are said to be blended.
        Units are in ARCSECONDS.

    n_proc : int, optional
        Number of processors to use. Should not exceed the number of cores.
        Default is one processor (no parallelization).

    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.
        
    hdf5_file_comp: str
        String of hdf5 file of companion events created in perform_pop_syn().
        Default is None.


    Returns
    -------
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
    _check_calc_events(hdf5_file, output_root2,
                       radius_cut, obs_time, n_obs, theta_frac,
                       blend_rad, n_proc, overwrite, hdf5_file_comp)

    ##########
    # Start of code
    #########

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
    
    if hdf5_file_comp is not None:
        hdc = itertools.repeat(hdf5_file_comp, reps)
        inputs = zip(llbb, hd, ot, no, rc, tf, br, hdc)
    
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

    if len(results_ev) > 0:
        events_tmp = np.concatenate(results_ev, axis=0)
        if len(results_bl) == 0:
            blends_tmp = np.array([])
        else:
            blends_tmp = np.concatenate(results_bl, axis=0)
        # Convert the events numpy recarray into an
        # Astropy Table for easier consumption.
        events_tmp = unique_events(events_tmp)
        events_final = Table(events_tmp)
        N_events = len(events_final)
        print('Candidate events detected: ', N_events)

        if len(results_bl) != 0:
            blends_tmp = unique_blends(blends_tmp)
        blends_final = Table(blends_tmp)
        
        #if len(results_cp) != 0:
         #   companions_tmp = unique_companions(companions_tmp, events_tmp)

        # Save out file
        events_final.write(output_root2 + '_events.fits', overwrite=overwrite)
        blends_final.write(output_root2 + '_blends.fits', overwrite=overwrite)
    else:
        N_events = 0
        print('No events!')
    t1 = time.time()

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    radius_cut = radius_cut / 1000.0  # back to arcsec
    popsycle_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                            cwd=popsycle_path).decode('ascii').strip()
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
    if hdf5_file_comp is not None:
        line8 += 'hdf5_file_comp , %s \n' % hdf5_file_comp

    line9 = 'VERSION INFORMATION' + '\n'
    line10 = str(now) + ' : creation date' + '\n'
    line11 = popsycle_hash + ' : PopSyCLE commit' + '\n'

    line12 = 'OTHER INFORMATION' + '\n'
    line13 = str(t1 - t0) + ' : total runtime (s)' + '\n'
    line14 = str(N_events) + ' : total number of events' + '\n'

    if N_events > 0:
        line15 = 'FILES CREATED' + '\n'
        line16 = output_root2 + '_events.fits : events file' + '\n'
        line17 = output_root2 + '_blends.fits : blends file' + '\n'
    else:
        line15 = 'NO FILES CREATED' + '\n'
        line16 = '\n'
        line17 = '\n'

    with open(output_root2 + '_calc_events.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3,
                        line4, line5, line6, line7, line8, empty_line,
                        line9, dash_line, line10, line11, empty_line,
                        line12, dash_line, line13, line14, empty_line, line15,
                        dash_line, line16, line17])

    print('calc_events runtime : {0:f} s'.format(t1 - t0))

    return


def _calc_event_time_loop(llbb, hdf5_file, obs_time, n_obs, radius_cut,
                          theta_frac, blend_rad, hdf5_file_comp = None):
    """
    Parameters
    -----------
    llbb : (int, int)
        Indices of (l,b) bin.

    obs_time, n_obs, radius_cut, theta_frac, blend_rad
    are all parameters of calc_events()

    Returns
    -------
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

    # Initialize events_llbb and blends_llbb and companions.
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

    # Adds separation in mas between primary and furthest companion if there are companions
    if hdf5_file_comp is not None:
        hfc = h5py.File(hdf5_file_comp, 'r')
        bigpatch_comp = np.hstack((hfc[name00], hfc[name01], hfc[name10], hfc[name11]))
        hfc.close()
        
        if len(bigpatch_comp) > 0:
            bigpatch_comp = rfn.append_fields(bigpatch_comp, 'sep', np.zeros(len(bigpatch_comp)), usemask = False) #separation in mas
            bigpatch_comp_df = pd.DataFrame(data = bigpatch_comp, columns = np.dtype(bigpatch_comp[0]).names)
            bigpatch_df = pd.DataFrame(data = bigpatch, columns = np.dtype(bigpatch[0]).names)

            rad = np.array(np.repeat(bigpatch_df['rad'], bigpatch_df['N_companions']))

            a_kpc = ((10**bigpatch_comp['log_a'])*unit.AU).to('kpc').value

            # abs(acos(i))
            bigpatch_comp_df['sep'] = np.abs(np.cos(bigpatch_comp['i']))*(np.arcsin(a_kpc/rad)*unit.radian).to('mas').value

            # Find max sep for triples
            bigpatch_comp_df_max = bigpatch_comp_df.groupby('system_idx').max()['sep'].reset_index()

            # Add separation to bigpatch
            sep = bigpatch_comp_df_max[['system_idx', 'sep']]
            # Cross referencing between separation from companions table and the primaries
            sep = sep.set_index("system_idx")
            bigpatch_df = bigpatch_df.set_index('obj_id')
            bigpatch_df = bigpatch_df.join(sep)
            bigpatch_df['sep'] = bigpatch_df['sep'].fillna(0) # Make nans from lack of companions to zeros
            
            bigpatch_df = bigpatch_df.reset_index() # Makes it index normally instead of by 'obj_id'
            bigpatch = rfn.append_fields(bigpatch, 'sep', bigpatch_df['sep'], usemask = False) 
            
            del bigpatch_df
            del bigpatch_comp_df
         
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
        theta_E = einstein_radius(bigpatch['systemMass'][lens_id],
                                  r_t[lens_id], r_t[sorc_id])  # mas      
        u = sep[event_id1] / theta_E
        
        binary_sep = None
        if 'sep' in bigpatch[0].dtype.names:
            binary_sep = bigpatch['sep'][event_id1]

        # Trim down to those microlensing events that really get close enough
        # to hope that we can detect them. Trim on a Theta_E criteria.
        event_lbt = _calc_event_cands_thetaE(bigpatch, theta_E, u, theta_frac,
                                             lens_id, sorc_id, time_array[i], binary_sep)

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
            # Note 3: Assumes all binaries are blended
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
    -----------
    bigpatch : array
        Compilation of 4 .h5 datasets containing stars.

    timei : float
        Time at which to evaluate.

    radius_cut : float
        Parameter of calc_events().
        Converted to mas

    Returns
    -------
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
                             sorc_id, timei, binary_sep = None):
    """
    Get sources and lenses that pass the radius cut.

    Parameters
    -----------
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

    Returns
    -------
    event_lbt : array
        Lenses and sources at a particular time t.

    """
    # If there are binaries extend the search radius to theta_frac + separation between primary
    # and furthest companion.
    if binary_sep is not None:
        theta_frac_comp = theta_frac + binary_sep
        if np.shape(u) != np.shape(theta_frac_comp):
            print(u, theta_frac_comp)
        #print(np.shape(u), np.shape(theta_frac_comp), np.shape(theta_frac), np.shape(bigpatch['sep']))
        adx = np.where(u < theta_frac_comp)[0]
    else: 
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
                                      u.value, usemask=False)
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
    -----------
    bigpatch : array
        Compilation of 4 .h5 datasets containing stars.

    c : SkyCoord object
        Coordinates of all the stars.

    event_lbt : array
        Lenses and sources at a particular time t.

    blend_rad : float
        Parameter of calc_events().

    Returns
    -------
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
    ----------
    event_table : numpy array 
        A table with all the events. There are equal numbers of columns
        containing info about the source and the lens, and four additional
        columns with theta_E, u0, mu_rel, and t0. The number of rows
        corresponds to the number of events.

    Returns
    -------
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
        min_idx = np.argmin(np.abs(dup_events['u0']))
        new_event_table[dpdx[ii]] = event_table[dup_idx[min_idx]]

    return new_event_table


def unique_blends(blend_table):
    """
    Given a blends table, there might be a single lens-source-neighbor triple
    multiple times (it might have been observed at different timesteps.)
    This function will eliminate duplicates, only keeping an event once. It
    is picked to be the first occurence.

    Parameters
    ----------
    blend_table : blend array 
        A table with all the events. There is 1 column with the unique
        source ID, 1 with the unique lens ID lens, 1 with the lens-neighbor
        separation, and the remaining columns contain info about the neighbors.


    Returns
    -------
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


def reduce_blend_rad(blend_tab, new_blend_rad, output_root, input_root = 'default', overwrite=False):
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
        
    input_root : str, optional
        The input root of the old blend table.
        If 'default', takes root of blend_tab.
    
    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exits the
        function if output files are already on disk.
        Default is False.
        

    Returns
    -------
    new_blend : .fits table
        New table with smaller blend radius.

    """
    if input_root == 'default':
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
    
    # Makes simlinks for necessary unchanged files
    links = ['_galaxia_params.txt', '_perform_pop_syn.log']
    for link in links:
        try:
            os.symlink(input_root + link, output_root + link)
        except:
            continue

    return




############################################################################
######### Refined event rate calculation and associated functions ##########
############################################################################


def _convert_photometric_99_to_nan(table, photometric_system='ubv'):
    for name in table.colnames:
        if ('exbv' in name) or (photometric_system in name):
            cond = np.where(table[name] == -99)[0]
            table[name][cond] = np.nan


def _check_refine_events(input_root, filter_name,
                         photometric_system, red_law, overwrite,
                         output_file, hdf5_file_comp, legacy, seed):
    """
    Checks that the inputs of refine_events are valid

    Parameters
    -----------
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
        The name of the reddening law to use from SPISEA.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        
    hdf5_file_comp: str
        String of hdf5 file of companion events created in perform_pop_syn().
        
    legacy : bool
        For running on files created before ~2020 when the filter system was changed
        to uppercase (i.e. from ubv_r to ubv_R) and before seeds were introduced.

    seed : None or int
        If not None, this forces the random orbit time for binaries to be fixed every time.
        If seed is added but there are no binaries, the seed will have no consequence.
    """

    if not isinstance(input_root, str):
        raise Exception('input_root (%s) must be a string.' % str(input_root))

    if not isinstance(filter_name, str):
        raise Exception('filter_name (%s) must be a string.' % str(filter_name))

    if not isinstance(photometric_system, str):
        raise Exception('photometric_system (%s) must be a string.' % str(photometric_system))

    if not isinstance(red_law, str):
        raise Exception('red_law (%s) must be a string.' % str(red_law))

    if not isinstance(output_file, str):
        raise Exception('output_file (%s) must be a string.' % str(output_file))

    if not isinstance(overwrite, bool):
        raise Exception('overwrite (%s) must be a boolean.' % str(overwrite))
        
    if not isinstance(hdf5_file_comp, str):
        if not isinstance(hdf5_file_comp, type(None)):
            raise Exception('hdf5_file_comp (%s) must be a str or a NoneType.' % str(hdf5_file_comp))
            
    if not isinstance(legacy, bool):
        raise Exception('legacy (%s) must be a boolean.' % str(legacy))

    if seed is not None:
        if not isinstance(seed, int):
            raise Exception('seed (%s) must be None or an integer.' % str(seed))

    # Check to see that the filter name, photometric system, red_law are valid
    if photometric_system not in photometric_system_dict:
        exception_str = 'photometric_system must be a key in ' \
                        'photometric_system_dict. \n' \
                        'Acceptable values are : '
        for photometric_system in photometric_system_dict:
            exception_str += '%s, ' % photometric_system
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    if filter_name not in photometric_system_dict[photometric_system]:
        exception_str = 'filter_name must be a value in ' \
                        'photometric_system_dict[%s]. \n' \
                        'Acceptable values are : ' % photometric_system
        for filter_name in photometric_system_dict[photometric_system]:
            exception_str += '%s, ' % filter_name
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    key = photometric_system + '_' + filter_name
    if red_law not in filt_dict[key]:
        exception_str = 'red_law must be a value in ' \
                        'filt_dict[%s]. \n' \
                        'Acceptable values are : ' % key
        for red_law in filt_dict[key]:
            exception_str += '%s, ' % red_law
        exception_str = exception_str[:-2]
        raise Exception(exception_str)


def refine_events(input_root, filter_name, photometric_system, red_law,
                  overwrite=False,
                  output_file='default', hdf5_file_comp=None, legacy = False, seed = None):
    """
    Takes the output Astropy table from calc_events, and from that
    calculates the time of closest approach. Will also return source-lens
    separation at this time.

    Parameters
    ----------
    input_root : str
        The root path and name of the \*_events.fits and \*_blends.fits.
        Don't include those suffixes yet.

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    photometric_system : str
        The name of the photometric system in which the filter exists.

    red_law : str
        The name of the reddening law to use from SPISEA.

    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    output_file : str, optional
        The name of the final refined_events file.
        If set to 'default', the format will be:
        <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits
            
    hdf5_file_comp: str or None, optional
        String of hdf5 file of companion events created in perform_pop_syn().
        Default is None.
    
    legacy : bool, optional
        For running on files created before ~2020 when the filter system was changed
        to uppercase (i.e. from ubv_r to ubv_R) and before seeds were introduced.
        Default is False.
        
    seed : int or None, optional
        If not None, this forces the random orbit time for binaries to be fixed every time.
        If seed is added but there are no binaries, the seed will have no consequence.
        Default is None.

    Returns
    -------
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits : fits file
        A file will be created that contains all the same objects, 
        only now with lots of extra columns of data.
    
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>_companions.fits : fits file
        If hdf5_file_comp is not None then a file will be created that contains the same objects 
        as a the companion file with extra columns of data.

    """
    # Set random seed
    np.random.seed(seed)
                      
    # Check if .fits file exists already. If it does, throw an error message
    # to complain and exit.
    if not overwrite and os.path.isfile(output_file):
        raise Exception('That refined_events.fits file name is taken! '
                        'Either delete the .fits file, or pick a new name.')

    # Error handling/complaining if input types are not right.
    _check_refine_events(input_root, filter_name,
                         photometric_system, red_law,
                         overwrite, output_file, hdf5_file_comp, legacy, seed)

    if output_file == 'default':
        output_file = '{0:s}_refined_events_{1:s}_{2:s}_{3:s}.fits'.format(input_root,
                                                                           photometric_system,
                                                                           filter_name,
                                                                           red_law)

    t_0 = time.time()

    event_fits_file = input_root + '_events.fits'
    blend_fits_file = input_root + '_blends.fits'
    galaxia_params_file = input_root + '_galaxia_params.txt'
    calc_events_log_file = input_root + '_calc_events.log'
    perform_pop_syn_log_file = input_root + '_perform_pop_syn.log'

    for filename in [event_fits_file,
                     blend_fits_file,
                     galaxia_params_file,
                     calc_events_log_file,
                     perform_pop_syn_log_file]:
        if not os.path.exists(filename):
            raise Exception(f'{filename} cannot be found.')

    event_tab = Table.read(event_fits_file)
    blend_tab = Table.read(blend_fits_file)
    
    # In legacy files set filters to upper case (i.e. ubv_r to ubv_R)
    if legacy == True:
        lowercase_ubv_filter_list = ['u', 'b', 'v', 'r', 'i', 'j', 'k', 'h']
        for ubv_filter in lowercase_ubv_filter_list:
            event_tab.rename_column('ubv_{}_S'.format(ubv_filter), 'ubv_{}_S'.format(ubv_filter.upper()))
            event_tab.rename_column('ubv_{}_L'.format(ubv_filter), 'ubv_{}_L'.format(ubv_filter.upper()))
            blend_tab.rename_column('ubv_{}_N'.format(ubv_filter), 'ubv_{}_N'.format(ubv_filter.upper()))
            
        

    # If photometric fields contain -99, convert to nan
    _convert_photometric_99_to_nan(event_tab, photometric_system)
    _convert_photometric_99_to_nan(blend_tab, photometric_system)

    # Only keep events with luminous sources
    event_tab = event_tab[~np.isnan(event_tab[photometric_system + '_' + filter_name + '_S'])]

    # Grab the obs_time from the calc_events log
    with open(calc_events_log_file, 'r') as my_file:
        for num, line in enumerate(my_file):
            if 'obs_time' in line.split(',')[0]:
                obs_time = line.split(',')[1]
                obs_time = float(obs_time)
                break

    # Grab the random seed from the galaxia param file
    with open(galaxia_params_file, 'r') as my_file:
        for num, line in enumerate(my_file):
            if 'seed' == line.split(' ')[0]:
                gal_seed = line.split(' ')[1].replace('\n', '')
                gal_seed = int(gal_seed)
                break

    # Grab the random seed from the perform_pop_syn log
    with open(perform_pop_syn_log_file, 'r') as my_file:
        for num, line in enumerate(my_file):

            if 'seed ' == line.split(',')[0]:
                pps_seed = line.split(',')[1].replace('\n', '')
                try:
                    pps_seed = int(pps_seed)
                except:
                    pps_seed = np.nan
                break
    
    # Sets random seed to nan for legacy files unless
    #  some value was set manually
    if legacy == True:
        try:
            gal_seed
        except NameError:
            gal_seed = np.nan
        
        try:
            pps_seed
        except NameError:
            pps_seed = np.nan
            

    # Calculate time and separation at closest approach, add to table
    # NOTE: calc_closest_approach modifies the event table!
    # It trims out events that peak outside the survey range!
    N_events_original = len(event_tab)
    print('Original candidate events: ', N_events_original)
    u0, t0 = calc_closest_approach(event_tab, obs_time)
    N_events_survey = len(event_tab)
    print('Candidate events in survey window: ', N_events_survey)
    event_tab['t0'] = t0  # days
    event_tab['u0'] = u0
    if len(event_tab) == 0:
        print('No events!')
        output_file = 'NO FILE CREATED'
    else:
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

        # Add random seeds as columns
        event_tab['pps_seed'] = np.ones(len(event_tab)) * pps_seed
        event_tab['gal_seed'] = np.ones(len(event_tab)) * gal_seed

        event_tab.write(output_file, overwrite=overwrite)
    
    

    
    
    # If Multiples
    if hdf5_file_comp is not None:
        hf_comp = h5py.File(hdf5_file_comp, 'r')
        companion_table = []
        
        #sets up event table to be indexed
        event_table = event_tab.to_pandas()
        event_tab_L = copy.deepcopy(event_table).set_index('obj_id_L')
        event_tab_S = copy.deepcopy(event_table).set_index('obj_id_S')
        # Loops through the non-metadata parts of hf_comp square by square
        # so the memory isn't overwhelmed
        for ii in np.hstack((hf_comp)):
            if ii != 'galaxyModelFile' and ii != 'lat_bin_edges' and ii != 'long_bin_edges':
                
                print('On square {}'.format(ii))
                
                patch_comp = np.hstack((hf_comp[ii],))
                if len(patch_comp) == 0:
                    continue
                    
                patch_comp_df = pd.DataFrame(data = patch_comp, columns = np.dtype(patch_comp[0]).names)
                
                # deletes patch so memory not overwhelmed
                del patch_comp
                
                #set the table to be indexed on system_idx
                #so that they can be joined based on those indices
                #event tables set to be indexed on obj_id_L, and obj_id_S above loop
                patch_comp_df = patch_comp_df.set_index('system_idx')
                
                patch_comp_df_columns = ['system_idx'] + list(patch_comp_df.columns) + ['obj_id_L', 'obj_id_S', 'prim_type']
                
                # only keeps columns where system is associated with a lens or source respectively
                patch_comp_L = patch_comp_df.join(event_tab_L, lsuffix='_comp', rsuffix='_prim', how='inner')
                patch_comp_S = patch_comp_df.join(event_tab_S, lsuffix='_comp', rsuffix='_prim', how='inner')
                
                # deletes patch so memory not overwhelmed
                del patch_comp_df
                
                # FIX ME: Memory could potentially be cleaned up here.
                # Reset index adds obj_id_L back as its own column
                # indexing on patch_comp_df_columns means that it only takes columns that were in the primary event table
                patch_comp_L['prim_type'] = "L"
                patch_comp_L['obj_id_L'] = (patch_comp_L.index).astype('float64')
                patch_comp_L['obj_id_S'] = (patch_comp_L['obj_id_S']).astype('float64')
                # In newer pandas versions name of index is maintained
                if patch_comp_L.index.name != 'system_idx':
                    patch_comp_L['system_idx'] = patch_comp_L.index
                if len(patch_comp_L) > 0:
                    patch_comp_L = patch_comp_L.reset_index()[patch_comp_df_columns]
                
                patch_comp_S['prim_type'] = "S"
                patch_comp_S['obj_id_S'] = (patch_comp_S.index).astype('float64')
                patch_comp_S['obj_id_L'] = (patch_comp_S['obj_id_L']).astype('float64')
                # In newer pandas versions name of index is maintained
                if patch_comp_S.index.name != 'system_idx':
                    patch_comp_S['system_idx'] = patch_comp_S.index
                if len(patch_comp_S) > 0:
                    patch_comp_S = patch_comp_S.reset_index()[patch_comp_df_columns]

                #For first time companion table is being created
                if len(companion_table) == 0:
                    if len(patch_comp_L) > 0:
                        companion_table = patch_comp_L
                        if len(patch_comp_S) > 0:
                            companion_table = pd.concat([companion_table, patch_comp_S], ignore_index=True)
                            #continue
                    elif len(patch_comp_S) > 0:
                        companion_table = patch_comp_S
                        #continue
                # For not first-time, appending       
                else:
                    # Appends them to companion table
                    # if-statement necessary because dtype of prim_type column doesn't change
                    # if the table is empty
                    if len(patch_comp_L) > 0:
                        companion_table = pd.concat([companion_table, patch_comp_L], ignore_index=True)
                    if len(patch_comp_S) > 0:
                        companion_table = pd.concat([companion_table, patch_comp_S], ignore_index=True)
                 
                # deletes patch so memory not overwhelmed
                del patch_comp_L
                del patch_comp_S
        hf_comp.close()
        
        
        if len(companion_table) > 0:
            
            # Adds parameters
            companion_table = _add_multiples_parameters(companion_table, event_tab)
            companion_table = _add_binary_angles(companion_table, event_tab)
            companion_table = Table.from_pandas(companion_table).filled(np.nan)
            
            # Adds metadata
            companion_table['i'].description = 'w/rt galactic galactic north'
            companion_table['Omega'].description = 'w/rt galactic galactic north'
            companion_table['omega'].description = 'w/rt galactic galactic north'
        
            # Writes fits file
            companion_table.write(output_file[:-5] + "_companions.fits", overwrite=overwrite)
            
    t_1 = time.time()

    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    popsycle_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=popsycle_path).decode('ascii').strip()
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
    line6 = popstar_hash + ' : SPISEA commit' + '\n'
    line7 = popsycle_hash + ' : PopSyCLE commit' + '\n'

    line8 = 'OTHER INFORMATION' + '\n'
    line9 = str(t_1 - t_0) + ' : total runtime (s)' + '\n'
    line10 = str(N_events_original) + ' : original candidate events (dark sources removed)' + '\n'
    line11 = str(N_events_survey) + ' : candidate events in survey window' + '\n'

    line12 = 'FILES CREATED' + '\n'
    line13 = output_file + ' : refined events'
    line14 = '\n' #By default no companion file created
    
    if hdf5_file_comp is not None:
        if len(companion_table) > 0:
            line14 = output_file[:-5] + "_companions.fits" + ' : companions refined events'

    with open(input_root + '_refined_events_' + photometric_system + '_' + filter_name + '_' + red_law + '.log', 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, empty_line,
                        line4, dash_line, line5, line6, line7, empty_line,
                        line8, dash_line, line9, line10, line11, empty_line,
                        line12, dash_line, line13, line14])

    print('refine_events runtime : {0:f} s'.format(t_1 - t_0))
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

    Returns
    -------
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
    With sign convention from Gould 2004:
    u0 > 0 then the source is to the east of the lens
    u0 < 0 then the source is to the west of the lens

    Parameters
    ----------
    event_tab : astropy table

    time : int/float

    Returns
    -------
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
     
    # Comment on sign conventions:
    # if u0_E > 0 then the source is to the East of the lens
    # if u0_E < 0 then the source is to the West of the lens
    # We adopt the following sign convention (same as Gould:2004):
    #    u0_amp > 0 means u0_E > 0
    #    u0_amp < 0 means u0_E < 0
    # Note that we assume beta = u0_amp (with same signs).

    sign = np.sign(c_S.icrs.ra - c_L.icrs.ra)
    
    return sign*u


def calc_blend_and_centroid(filter_name, red_law, blend_tab, photometric_system='ubv'):
    """
    Given the absolute magnitudes of a bunch of things,
    calculate their blended apparent magnitude and flux.
    Also calculate the centroid of these things.
    Assumes all binaries are blended.

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
    if type(flux_N) == np.ma.core.MaskedArray or type(flux_N) == MaskedColumn:
        flux_N = np.nan_to_num(flux_N.filled(np.nan))
    else:
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
    -----------
    filter_name : str

    event_tab : Astropy table

    blend_tab : Astropy table

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
    
    if type(flux_L) == np.ma.core.MaskedArray or type(flux_L) == MaskedColumn:
        flux_L = np.nan_to_num(flux_L.filled(np.nan))
    else:
        flux_L = np.nan_to_num(flux_L)

    if type(flux_S) == np.ma.core.MaskedArray or type(flux_S) == MaskedColumn:
        flux_S = np.nan_to_num(flux_S.filled(np.nan))
    else:
        flux_S = np.nan_to_num(flux_S)

    ##########
    # Find the blends.
    ##########
    # LS_pairs are the (unique) array of lens/source pairs that correspond
    #     to a single microlensing event.
    # blend_pairs is a list of all lens/source pairs associated with a 
    #     neighbor. A lens/source pair is listed multiple times, listed N
    #     times if there are N neighbors.
    LS_pairs = np.stack((event_tab['obj_id_L'].data, event_tab['obj_id_S'].data),
                        axis=-1)
    blend_pairs = np.stack((blend_tab['obj_id_L'].data, blend_tab['obj_id_S'].data),
                           axis=-1)

    # flux_N will store the total neighbor flux associated with an event
    #     lens/source pair.
    # cent_glon_N and cent_glat_N will store the centroid.
    flux_N = np.zeros(len(app_L))
    event_tab['cent_glon_' + filter_name + '_N'] = np.zeros(len(app_L))
    event_tab['cent_glat_' + filter_name + '_N'] = np.zeros(len(app_L))

    # Get the unique blend_pairs and the first instance of each.
    # The blend pairs are listed sequentially. That is, if you have three
    #     events E1, E2, E3 which have 3, 5, and 2 neighbors respectively,
    #     they will be listed in the table as 
    #         N1_1, N1_2, N1_3, N2_1, N2_2, N2_3, N2_4, N2_5, N3_1, N3_2
    #     and the corresponding unique indices uni_bidx will be
    #         0, 3, 8.
    #     We then append the length of the blend_pairs table (in this case 10)
    #     to uni_bidx to get uni_bidxx. That way, you can get all the indices
    #     of a particular event by indexing from
    #     uni_bidxx[idx][0] to uni_bidxx[idx+1][0].
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
    delta_m = calc_bump_amp(np.abs(event_tab['u0']), flux_S, flux_L, flux_N)
    event_tab['delta_m_' + filter_name] = delta_m

    # Calculate the blend fraction
    # (commonly used in microlensing): f_source / f_total
    f_blend = flux_S / flux_tot
    event_tab['f_blend_' + filter_name] = f_blend

    return


def _add_binary_angles(companion_table, event_table):
    """
    -Adds the following angles to the companion array:
        - alpha: angle between North and binary axis (East of North)
        - phi_pi_E: angle between North and proper motion vector (East of North)
        - phi: angle between the proper motion and the binary axis
    -Each set of angles is for the source if the companion is 
    assocaited with the source or the lens if the companion is associated with the lens.
    -Binary parameters established in a galactic spherical coordinate systems.
    -Big Omega measured from galactic north increasing in the direction of galactic east (positive l).
    
    Parameters
    -----------
    companion_table : pandas table
        Pandas table from companion events which were matched from the events table.
        
    event_table : astropy table
        Astropy table with refined events

    Returns
    -------
    companion_tmp_df : pandas table
        Companions table modified to include alpha, phi_pi_E, phi.
    
    
    """
    companion_tmp_df = copy.deepcopy(companion_table)    
    event_table_df = event_table.to_pandas()
    
    #indexes on both obj_id_L and obj_id_S at once
    companion_tmp_df = companion_tmp_df.set_index(['obj_id_L', 'obj_id_S'])
    event_table_df = event_table_df.set_index(['obj_id_L', 'obj_id_S'])
    companion_tmp_df_joined = companion_tmp_df.join(event_table_df, lsuffix='_comp', rsuffix='_prim', how='inner')
        
    del event_table_df
    
    start_time_bin_angles = time.time()    
    alphas, phi_pi_Es, phis = calculate_binary_angles(companion_tmp_df_joined)
        
    companion_tmp_df_joined['alpha'] = alphas
    companion_tmp_df_joined['phi_pi_E'] = phi_pi_Es
    companion_tmp_df_joined['phi'] = phis
    
    # reset index adds obj_id_L and obj_id_S back as their own columns
    # doing list(cluster.companions.columns) + those other names means that it only takes columns 
    # that were in cluster.companions and obj_id_L, obj_id_S, alpha,...
    companion_tmp_df = companion_tmp_df_joined.reset_index()[list(companion_tmp_df.columns) + ['obj_id_L', 'obj_id_S'] + ['alpha', 'phi_pi_E', 'phi']]
    #print('full bin angles loop', time.time() - start_time_bin_angles)
    
    del companion_tmp_df_joined
    return companion_tmp_df

def calculate_binary_angles(joined_table):
    """
    * Calculates the following angles to the companion array:
    
        * alpha: angle between North and binary axis (East of North)
        * phi_pi_E: angle between North and proper motion vector (East of North)
        * phi: angle between the proper motion and the binary axis
    
    * Each set of angles is for the source if the companion is 
      assocaited with the source or the lens if the companion is associated with the lens.
    * Binary parameters established in a galactic spherical coordinate systems.
    * Big Omega measured from galactic north increasing in the direction of galactic east (positive l).
    
    Parameters
    ----------
    joined_table : pandas table
        Joined companion and event table on obj_id_L and obj_id_S

    Returns
    -------
    alphas : list
        Angle between North and binary axis (East of North)
    phi_pi_Es : list
        Angle between North and proper motion vector (East of North)
    phis : list
        Angle between the proper motion and the binary axis (phi_pi_E - alpha)
    
    
    """
    alphas = []
    phi_pi_Es = []
    phis = []
    
    for index, row in joined_table.iterrows():
        
        # At the end of refine events prim_type is encoded in bytes
        # so we need to decode it to strings
        # if we feed in a table from and alternate location
        prim_type = row['prim_type']
        if type(prim_type) == bytes:
            prim_type = prim_type.decode("utf-8") 
            
        if prim_type == 'L':
            not_prim_type = 'S'
        else:
            not_prim_type = 'L'
        
        # First calculate binary axis
        orb = orbits.Orbit()
        orb.w = row['omega'] # [degrees]
        orb.o = row['Omega'] # [degrees]
        orb.i = row['i'] # [degrees]
        orb.e = row['e'] # float between 0 and 1
        orb.p = row['P'] # [years]
        orb.t0 = np.random.rand()*row['P'] + row['t0'] # [years] This is initial
        orb.mass = row['systemMass_{}'.format(prim_type)] # [Msun]
        
        
        # Position of the companion when primary at origin
        # r[0][0] is galactic east
        # r[0][1] is galactic north
        # r[0][2] is line of sight away from the sun
        # Assumes the primary is at the origin (0,0)
        
        (r, v, a) = orb.kep2xyz(np.array([row['t0']]))

        r_kpc = np.array([(ii*unit.AU).to('kpc').value for ii in r[0]])
        
        r_mas = (np.arcsin(r_kpc[0:2]/row['rad_{}'.format(prim_type)])*unit.radian).to('mas').value
        
        # Change in positional offset (b and lcosb) from companion to the primary 
        # pointing to companion (assuming primary at 0,0)
        deltaz_lcosb = r_mas[0]
        deltaz_b = r_mas[1]
        
        #b and l of primary
        z1_b = row['glat_{}'.format(prim_type)]*unit.degree
        z1_l = row['glon_{}'.format(prim_type)]*unit.degree
        
        
        deltaz_l = deltaz_lcosb/np.cos(np.deg2rad(z1_b)) #Need to find the spherical l instead of projected for SkyCoord
        
        #add positions of the primary to the deltas to find b and l of companion
        z2_b = (deltaz_b*unit.mas + z1_b.to('mas')).to('degree')
        z2_l = (deltaz_l*unit.mas + z1_l.to('mas')).to('degree')
        
        #b and l of primary of system companion not part of
        not_prim_b = row['glat_{}'.format(not_prim_type)]*unit.degree
        not_prim_l = row['glon_{}'.format(not_prim_type)]*unit.degree
        
        #SkyCoord objects for the primary, companion, and primary of system companion not part of
        coord_object_z1 = SkyCoord(l = z1_l, b = z1_b, 
                                   pm_l_cosb = row['mu_lcosb_{}'.format(prim_type)]*unit.mas/unit.year, 
                                   pm_b = row['mu_b_{}'.format(prim_type)]*unit.mas/unit.year, 
                                   frame ='galactic')
        
        coord_object_z2 = SkyCoord(l = z2_l, b = z2_b, 
                                   pm_l_cosb = row['mu_lcosb_{}'.format(prim_type)]*unit.mas/unit.year, 
                                   pm_b = row['mu_b_{}'.format(prim_type)]*unit.mas/unit.year, 
                                   frame ='galactic')
        
        coord_object_not_prim = SkyCoord(l = not_prim_l, b = not_prim_b, 
                                         pm_l_cosb = row['mu_lcosb_{}'.format(not_prim_type)]*unit.mas/unit.year, 
                                         pm_b = row['mu_b_{}'.format(not_prim_type)]*unit.mas/unit.year, frame ='galactic')

        #ra and dec of primary, find racosdec to find projected coords
        z1_ra = coord_object_z1.icrs.ra.value
        z1_dec = coord_object_z1.icrs.dec.value
        z1_racosdec = z1_ra*np.cos(np.deg2rad(z1_dec))
        
        #ra and dec of companion, find racosdec to find projected coords
        z2_ra = coord_object_z2.icrs.ra.value
        z2_dec = coord_object_z2.icrs.dec.value
        z2_racosdec = z2_ra*np.cos(np.deg2rad(z2_dec))
        
        #Find change in dec and racosdec companion - primary (so vector is pointed at companion)
        deltaz_dec = z2_dec - z1_dec
        deltaz_racosdec = (z2_ra - z1_ra)*np.cos(np.deg2rad(z1_dec)) #z2_racosdec - z1_racosdec
        
        
        #Proper motions of the primary, comanion, and source in racosdec and dec
        mu_racosdec_z1 = coord_object_z1.icrs.pm_ra_cosdec.value
        mu_dec_z1 = coord_object_z1.icrs.pm_dec.value
        
        mu_racosdec_z2 = coord_object_z2.icrs.pm_ra_cosdec.value
        mu_dec_z2 = coord_object_z2.icrs.pm_dec.value
        
        mu_racosdec_not_prim = coord_object_not_prim.icrs.pm_ra_cosdec.value
        mu_dec_not_prim = coord_object_not_prim.icrs.pm_dec.value
        
        #Relative proper motion between the primary and companion
        #delta_mu_racosdec.append(np.abs(mu_racosdec_z1 - mu_racosdec_z2))
        #delta_mu_dec.append(np.abs(mu_dec_z1 - mu_dec_z2))
        
        if np.abs(mu_racosdec_z1 - mu_racosdec_z2) > 0.01 or np.abs(mu_dec_z1 - mu_dec_z2) > 0.01:
            print('***************** WARNING ******************')
            print('Discrepancy between companion and primary proper motion components > 0.01 mas/year')
            print('********************************************')
            print('mu_racosdec_z1', mu_racosdec_z1)
            print('mu_racosdec_z2', mu_racosdec_z2)
            print('mu_dec_z1', mu_dec_z1)
            print('mu_dec_z2', mu_dec_z2)

        #Relative proper motion between primary and system companion not part of
        mu_racosdec_rel = mu_racosdec_not_prim - mu_racosdec_z1
        mu_dec_rel = mu_dec_not_prim - mu_dec_z1
        
        
        # Since deltaz is pointed toward the companion, subtract 180 degrees mod 360 to find alpha
        alpha = (np.rad2deg(np.arctan2(deltaz_racosdec,deltaz_dec)) - 180) % 360
        
        # Mod 360 so it'll be in the same range as alpha
        phi_pi_E = np.rad2deg(np.arctan2(mu_racosdec_rel,mu_dec_rel)) % 360
        
        phi = (alpha - phi_pi_E) % 360
        
        alphas.append(alpha)
        phi_pi_Es.append(phi_pi_E)
        phis.append(phi)
    
    return alphas, phi_pi_Es, phis

def _add_multiples_parameters(companion_table, event_table):
    """
    Adds mass ratio, separation in mas between primary and companion,
    and period in years to companions array.

    Parameters
    -----------
    companion_table : pandas table
        Pandas table from companion events which were matched from the events table.
        
    event_table : astropy table
        Astropy table with refined events

    Returns
    -------
    companion_tmp_df : pandas table
        Companions table modified to include new parameters (q, sep, and P).

    """
    companion_tmp_df = copy.deepcopy(companion_table)
    event_table_df = event_table.to_pandas()
    
    
    #indexes on both obj_id_L and obj_id_S at once
    companion_tmp_df = companion_tmp_df.set_index(['obj_id_L', 'obj_id_S'])
    event_table_df = event_table_df.set_index(['obj_id_L', 'obj_id_S'])
    companion_tmp_df_joined = companion_tmp_df.join(event_table_df, lsuffix='_comp', rsuffix='_prim', how='inner')

    del event_table_df
    
    a_kpc = (np.array(10**companion_tmp_df_joined['log_a'])*unit.AU).to('kpc').value
    
    event_prim_masses = np.zeros(len(companion_tmp_df_joined))
    event_system_masses = np.zeros(len(companion_tmp_df_joined))
    event_prim_distances = np.zeros(len(companion_tmp_df_joined))
    counter = np.arange(len(companion_tmp_df_joined))
    for row, count in zip(companion_tmp_df_joined.iterrows(), counter):
        # indexing row[1] is because row[0] is the row's index
        prim_type = row[1]['prim_type']
        
        # this is activated only for testing, since prim_type is saved
        # as a byte in the final fits file
        if type(prim_type) == bytes:
            prim_type = prim_type.decode("utf-8") 
            
        event_prim_masses[count] = row[1]['mass_' + prim_type]
        event_system_masses[count] = row[1]['systemMass_' + prim_type]
        event_prim_distances[count] = row[1]['rad_' + prim_type]
        
    companion_tmp_df_joined['q'] = companion_tmp_df_joined['mass']/event_prim_masses #mass ratio
    #abs(acos(i))
    abs_cos_i = np.abs(np.cos(list(companion_tmp_df_joined['i'])))
    a_mas = (np.arcsin(a_kpc/event_prim_distances)*unit.radian).to('mas').value
    companion_tmp_df_joined['sep'] = abs_cos_i*a_mas #separation in mas
    companion_tmp_df_joined['P'] = orbits.a_to_P(event_system_masses, 10**companion_tmp_df_joined['log_a']) #period
    # reset index adds obj_id_L and obj_id_S back as their own columns
    # doing list(cluster.companions.columns) + those other names means that it only takes columns 
    # that were in cluster.companions and obj_id_L, obj_id_S, q,...
    companion_tmp_df = companion_tmp_df_joined.reset_index()[list(companion_tmp_df.columns) + ['obj_id_L', 'obj_id_S'] + ['q', 'sep', 'P']]

    del companion_tmp_df_joined
    return companion_tmp_df


############################################################################
##### Refined binary event rate calculation and associated functions #######
############################################################################


def _check_refine_binary_events(events, companions,
                                photometric_system, filter_name,
                                overwrite, output_file,
                                save_phot, phot_dir, n_proc, multi_proc):
    """
    Checks that the inputs of refine_binary_events are valid

    Parameters
    -----------
    events : str
        fits file containing the events calculated from refine_events
    
    companions : str
        fits file containing the companions calculated from refine_events
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.

    output_file : str
        The name of the final refined_events file.
        If set to 'default', the format will be:
            <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits
            
    save_phot : bool
        If set to True, saves the photometry generated instead of just parameters.
    
    phot_dir : str
        Name of the directory photometry is saved if save_phot = True.
        This parameters is NOT optional if save_phot = True.
        
    n_proc : int
        Number of processors to use. Should not exceed the number of cores.
        Default is one processor (no parallelization).
    
    multi_proc : bool
        Even if n_proc = 1, a pool is still created. If multi_proc = False,
        instead there is just a for-loop to generate and analyze the lightcurves.
        If multi_proc == False, n_proc must = 1.
    """

    if not isinstance(events, str):
        raise Exception('events (%s) must be a string.' % str(events))

    if not isinstance(companions, str):
        raise Exception('companions (%s) must be a string.' % str(companions))
    
    if not isinstance(photometric_system, str):
        raise Exception('photometric_system (%s) must be a string.' % str(photometric_system))
    
    if not isinstance(filter_name, str):
        raise Exception('filter_name (%s) must be a string.' % str(filter_name))

    if not isinstance(output_file, str):
        raise Exception('output_file (%s) must be a string.' % str(output_file))

    if not isinstance(overwrite, bool):
        raise Exception('overwrite (%s) must be a boolean.' % str(overwrite))
        
    if not isinstance(save_phot, bool):
        raise Exception('save_phot (%s) must be a boolean.' % str(save_phot))
        
    if not isinstance(phot_dir, str):
        if not isinstance(phot_dir, type(None)):
            raise Exception('phot_dir (%s) must be a string or None.' % str(phot_dir))
            
    if not isinstance(n_proc, int):
        raise Exception('n_proc (%s) must be an integer.' % str(seed))
    
    if not isinstance(multi_proc, bool):
        raise Exception('multi_proc (%s) must be a boolean.' % str(multi_proc))
    
    if n_proc > 1 and multi_proc == False:
        raise Exception('if multi_proc is False, n_proc must = 1')

    # Check to see that the filter name, photometric system, red_law are valid
    if photometric_system not in photometric_system_dict:
        exception_str = 'photometric_system must be a key in ' \
                        'photometric_system_dict. \n' \
                        'Acceptable values are : '
        for photometric_system in photometric_system_dict:
            exception_str += '%s, ' % photometric_system
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

    if filter_name not in photometric_system_dict[photometric_system]:
        exception_str = 'filter_name must be a value in ' \
                        'photometric_system_dict[%s]. \n' \
                        'Acceptable values are : ' % photometric_system
        for filter_name in photometric_system_dict[photometric_system]:
            exception_str += '%s, ' % filter_name
        exception_str = exception_str[:-2]
        raise Exception(exception_str)

def refine_binary_events(events, companions, photometric_system, filter_name,
                         red_law = 'Damineli16',
                         overwrite = False, output_file = 'default',
                         save_phot = False, phot_dir = None,
                         n_proc = 1, multi_proc = True):
    """
    Takes the output Astropy table from refine_events 
    (both primaries and companions) and from that calculates 
    the binary light curves.

    Parameters
    ----------
    events : str
        fits file containing the events calculated from refine_events
    
    companions : str
        fits file containing the companions calculated from refine_events
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    red_law : str, optional
        Reddening law. Default is Damineli16
        
    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    output_file : str, optional
        The name of the final refined_events file.
        If set to 'default', the format will be:
        <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits
            
    save_phot : bool, optional
        If set to True, saves the photometry generated instead of just parameters.
        Default is False.
    
    phot_dir : str or None, optional
        Name of the directory photometry is saved if save_phot = True.
        This parameters is NOT optional if save_phot = True.
        Default is None.
        
    n_proc : int, optional
        Number of processors to use. Should not exceed the number of cores.
        Default is one processor (no parallelization).
    
    multi_proc : bool, optional
        Even if n_proc = 1, a pool is still created. If multi_proc = False,
        instead there is just a for-loop to generate and analyze the lightcurves.
        If multi_proc == False, n_proc must = 1.
        Default is True.
    
    Returns
    -------
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>_companions_rb.fits : fits file
        A file that contains all the same objects, only now with lots of extra
        columns of data. (rb stands for refine binaries).
    
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>_companions_rb_mp.fits : fits file
        A file will be created named that contains the data for each individual 
        peak for events with multiple peaks.
        (mp stands for multiple peaks).

    """
    start_time = time.time()
    
    if not overwrite and os.path.isfile(output_file):
        raise Exception('That refined_events.fits file name is taken! '
                        'Either delete the .fits file, or pick a new name.')

    if save_phot == True and phot_dir == None:
        raise Exception('phot_dir is "none". Input a directory to save photometry.')
        
    # Error handling/complaining if input types are not right.
    _check_refine_binary_events(events, companions, 
                         photometric_system, filter_name,
                         overwrite, output_file,
                         save_phot, phot_dir, n_proc, multi_proc)
        

    event_table = Table.read(events)
    comp_table = Table.read(companions)

    comp_table['companion_idx'] = np.arange(len(comp_table))

    comp_table.rename_column('m_ukirt_H', 'm_ubv_H')
    comp_table.rename_column('m_ukirt_J', 'm_ubv_J')
    comp_table.rename_column('m_ukirt_K', 'm_ubv_K')

    event_table['f_blend_%s' % filter_name] = event_table['f_blend_%s' % filter_name] # None of these should be nan
    if type(comp_table['m_%s_%s' % (photometric_system, filter_name)]) == np.ma.core.MaskedArray or type(comp_table['m_%s_%s' % (photometric_system, filter_name)]) == MaskedColumn:
        comp_table['m_%s_%s' % (photometric_system, filter_name)] = comp_table['m_%s_%s' % (photometric_system, filter_name)].filled(np.nan)
    
    event_table.add_column( Column(np.zeros(len(event_table), dtype=float), name='n_peaks') )
    event_table.add_column( Column(np.zeros(len(event_table), dtype=float), name='bin_delta_m') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='tE_sys') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='tE_primary') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='primary_t') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='avg_t') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='std_t') )
    event_table.add_column( Column(np.empty(len(event_table), dtype=float), name='asymmetry') )
    
    event_table['bin_delta_m'][:] = np.nan
    event_table['tE_sys'][:] = np.nan
    event_table['tE_primary'][:] = np.nan
    event_table['primary_t'][:] = np.nan
    event_table['avg_t'][:] = np.nan
    event_table['std_t'][:] = np.nan
    event_table['asymmetry'][:] = np.nan
    
    lightcurve_table = Table(names=('obj_id_L', 'obj_id_S', 'companion_id_L', 'companion_id_S', 'class', 'n_peaks', 'bin_delta_m', 
                                     'tE_sys', 'tE_primary', 'primary_t', 'avg_t', 'std_t', 'asymmetry', 'mp_rows', 'used_lightcurve'),
                            dtype = (float, float, float, float, str, float, float, float, float, float, float, float, float, dict, bool))
    
    # This table is for events with more than one peak to characterize those peaks
    # comp_id is position of companion in companion table for reference
    # obj_id_L, obj_id_S, and n_peaks are the same as in the companion table, just for reference
    # t is time of peak
    # tE is Einstein crossing time of peak defined by times of 0.5*(max(peak mag) - min(peak mag))
    # delta m is the change in magnitude between the peak and baseline
    # ratio is the magnitude ratio between min peak/max peak
    mult_peaks = Table(names=('companion_idx', 'obj_id_L', 'obj_id_S', 'n_peaks', 't', 'tE', 'delta_m', 'ratio'),
                       dtype = (object, float, float, float, float, float, float, float))
    
    multiples_lightcurves = sum((event_table['isMultiple_S'] == 1) | (event_table['isMultiple_L'] == 1))
    
    # Set up the multiprocessing
    if multi_proc:
        pool = Pool(n_proc)
    
    grouped_comps = comp_table.group_by(['obj_id_L', 'obj_id_S'])
    
    event_table_df = event_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'])
    empty_lists = [ [] for _ in range(len(event_table_df)) ]
    event_table_df['companion_idx_list'] = empty_lists
    
    inputs = np.empty(multiples_lightcurves, dtype = object)
    for i in range(len(grouped_comps.groups)):
        obj_id_L = grouped_comps.groups.keys[i][0]
        obj_id_S = grouped_comps.groups.keys[i][1]
        obj_id_L_S = (obj_id_L, obj_id_S)
        event_table_df['companion_idx_list'].loc[obj_id_L_S] = list(grouped_comps.groups[i]['companion_idx'])
        inputs[i] = [[event_table_df.loc[obj_id_L_S]], grouped_comps.groups[i].to_pandas(), obj_id_L, obj_id_S, 
                     photometric_system, filter_name, red_law, save_phot, phot_dir, overwrite]
    
    one_lightcurve_analysis
    
    if multi_proc:
        results = pool.starmap(one_lightcurve_analysis, inputs)
        # Multi-Threaded: clean up pool
        pool.close()
        pool.join()
    elif multi_proc == False:
        results = np.zeros(len(inputs), dtype=object)
        for i, input in enumerate(inputs):
            results[i] = one_lightcurve_analysis(*input)
    
    result_keys = list(results[0][0][0].keys())
    for i in range(len(grouped_comps.groups)):
        obj_id_L = grouped_comps.groups.keys[i][0]
        obj_id_S = grouped_comps.groups.keys[i][1]
        obj_id_L_S = (obj_id_L, obj_id_S)
        for result in results[i]:
            if result[0]['used_lightcurve'] == 1:
                for key in result_keys:
                    # only lightcurves with highest delta m (if there were triples and
                    # multiple were simulated) go into event_table
                    if key != 'mp_rows' and key != 'used_lightcurve':
                        event_table_df[key].loc[obj_id_L_S] = result[0][key]
            
            try:
                lightcurve_table.add_row(result[1] | result[0])
            except TypeError:
                # Fails when using python version < 3.9
                lightcurve_table.add_row({**result[1], **result[0]})
            
            rows = result[0]['mp_rows']
            for row in rows:
                mult_peaks.add_row(row)
    
    event_table_df['companion_idx_list'] = event_table_df['companion_idx_list'].astype(str)
    event_table = Table.from_pandas(event_table_df.reset_index()).filled(np.nan)
    
    lightcurve_table.remove_column('mp_rows')

    # For write issues, if mult_peaks is empty, define it without dtype
    if len(mult_peaks) == 0:
        mult_peaks = Table(names=('companion_idx', 'obj_id_L', 'obj_id_S', 'n_peaks', 't', 'tE', 'delta_m', 'ratio'))


    # Writes fits file
    if output_file == 'default':
        event_table.write(events[:-5] + "_rb.fits", overwrite=overwrite)
        comp_table.write(companions[:-5] + "_rb.fits", overwrite=overwrite)
        mult_peaks.write(companions[:-5] + "_rb_multi_peaks.fits", overwrite=overwrite)
        lightcurve_table.write(events[:-5] + "_rb_lightcurves.fits", overwrite=overwrite)
    else:
        event_table.write(output_file + ".fits", overwrite=overwrite)
        comp_table.write(output_file + "_companions.fits", overwrite=overwrite)
        mult_peaks.write(output_file + "_multi_peaks.fits", overwrite=overwrite)
        lightcurve_table.write(output_file + "_lightcurves.fits", overwrite=overwrite)
        
        
        
    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    popsycle_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=popsycle_path).decode('ascii').strip()
    popstar_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                           cwd=popstar_path).decode('ascii').strip()
    
    end_time = time.time()
    
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'event_file : ' + events + '\n'
    line2 = 'companion_file : ' + companions + '\n'
    line3 = 'save_phot : ' + str(save_phot) + '\n'
    line3p5 = 'n_proc : ' + str(n_proc) + '\n'

    line4 = 'VERSION INFORMATION' + '\n'
    line5 = str(now) + ' : creation date' + '\n'
    line6 = popstar_hash + ' : SPISEA commit' + '\n'
    line7 = popsycle_hash + ' : PopSyCLE commit' + '\n'

    line8 = 'OTHER INFORMATION' + '\n'
    line9 = str(end_time - start_time) + ' : total runtime (s)' + '\n'
    line10 = str(len(lightcurve_table)) + ' : number of simulated lightcurves' + '\n'

    line11 = 'FILES CREATED' + '\n'
    if output_file == 'default':
        line11p1 = events[:-5] + "_rb.fits" + ' : event table w/ binary params' + '\n'
        line12 = companions[:-5] + "_rb.fits" + ' : companion table w/ idx' + '\n'
        line13 = companions[:-5] + "_rb_multi_peaks.fits" + ' : multiple peak events table' + '\n'
        line13p5 = events[:-5] + "_rb_lightcurves.fits" + ' : lightcurve file w/ included+non-included lightcurves' + '\n'
        log_name = companions[:-5] + "_rb.log"
    else:
        line11p1 = output_file + ".fits" + ' : event table w/ binary params' + '\n'
        line12 = output_file + "_companions.fits" + ' : companion table w/ idx' + '\n'
        line13 = output_file + "_multi_peaks.fits" + ' : multiple peak events table' + '\n'
        line13p5 = output_file + "_lightcurves.fits" + ' : lightcurve file w/ included+non-included lightcurves' + '\n'
        log_name = output_file + '.log'
     
    line14 = '\n'
    if save_phot == True:
        line14 = phot_dir + ' : directiory of photometry'
       
    
    with open(log_name, 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, line3p5, empty_line,
                        line4, dash_line, line5, line6, line7, empty_line,
                        line8, dash_line, line9, line10, empty_line,
                        line11, dash_line, line11p1, line12, line13, line13p5, line14])

    print('refine_binary_events runtime : {0:f} s'.format(end_time - start_time))
    return

def one_lightcurve_analysis(event_table_row, comp_table_rows, obj_id_L, obj_id_S,
                            photometric_system, filter_name, red_law,
                            save_phot = False, phot_dir = None, overwrite = False):
    """
    Generate BAGLE model, photometry, and generate binary lightcurve parameters.
    
    Parameters
    ----------
    event_table_row : Astropy table
        Astropy table row from event table with event parameters.
    
    comp_table_rows: Astropy table
        Astropy table with row(s) associated with event row. I.e. if it's a binary source,
        binary lens, there should be two rows.
    
    obj_id_L : int
        Object id of the lens associated with event.
    
    obj_id_S : int
        Object id of the source associated with event.
    
    photometric_system : str
        Name of the photometric system, i.e. 'ubv'.
    
    filter_name : str
        Name of filter associated with photometric system, i.e. 'I'.
    
    red_law : str
        Name of reddening law in filt_dict list above, i.e. 'Damineli16'.
    
    save_phot : bool, optional
        If set to True, saves the photometry generated instead of just parameters.
        Default is False.
    
    phot_dir : str or None, optional
        Name of the directory photometry is saved if save_phot = True.
        This parameters is NOT optional if save_phot = True.
        Default is None.
        
    overwrite : bool, optional
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.
    
    Returns
    -------
    lightcurve_parameters : list
        List of two dictionaries:
            - Dictionary of additional parameters added to the event table
              and one called 'mp_rows' which are the rows to be
              added to the multi peak table.
            - Dictonary of lightcurve parameters which include both lightcurves
              whose parameters are saved in the event table and those that are not.
              i.e. in a triple system only the one with the highest delta_m's p
              arameters are used for the lightcurve table. These will all be saved in
              a separate lightcurve_table.
    """
    
    multi_L = event_table_row[0]['isMultiple_L']
    multi_S = event_table_row[0]['isMultiple_S']
    if multi_L == 1 and multi_S == 1:
        event_type = 'BSBL'
    elif multi_L == 1 and multi_S == 0:
        event_type = 'PSBL'
    elif multi_L == 0 and multi_S == 1:
        event_type = 'BSPL'
    else:
        raise Exception('one_lightcurve_analysis() only analyizes binary events')

    lightcurve_parameters = []
    max_delta_m = np.nan
    if event_type == 'BSBL':
        comp_idxs_S = np.where(comp_table_rows['prim_type'] == b"S")[0]
        comp_idxs_L = np.where(comp_table_rows['prim_type'] == b"L")[0]
        for comp_idx_S in comp_idxs_S:
            global_comp_idx_S = comp_table_rows['companion_idx'][comp_idx_S]
            for comp_idx_L in comp_idxs_L:
                global_comp_idx_L = comp_table_rows['companion_idx'][comp_idx_L]
                name = "L_{}_S_{}".format(obj_id_L, obj_id_S) + "compL_{}_compS_{}".format(global_comp_idx_L, global_comp_idx_S)
                model_parameter_dict, _, _ = get_bsbl_lightcurve_parameters(event_table_row, comp_table_rows, int(comp_idx_L), int(comp_idx_S), 
                                                                      photometric_system, filter_name, red_law, event_id = 0)
                model = bsbl_model_gen(model_parameter_dict)
                param_dict = lightcurve_parameter_gen(model, model_parameter_dict, np.array([global_comp_idx_L, global_comp_idx_S]), 
                                                      obj_id_L, obj_id_S, name, save_phot, phot_dir, overwrite)
                lightcurve_dict = {'obj_id_L' : obj_id_L, 'obj_id_S' : obj_id_S, 'companion_id_L' : comp_table_rows['companion_idx'][comp_idx_L], 
                                   'companion_id_S' : comp_table_rows['companion_idx'][comp_idx_S], 'class' : event_type}
                lightcurve_parameters.append([param_dict, lightcurve_dict])
                
            del global_comp_idx_S, global_comp_idx_L, param_dict,
            lightcurve_dict, model, model_parameter_dict, name
            
        del comp_idxs_S, comp_idxs_L
                
    elif event_type == 'PSBL':
        for comp_idx in range(len(comp_table_rows)):
            global_comp_idx = comp_table_rows['companion_idx'][comp_idx]
            name = "L_{}_S_{}".format(obj_id_L, obj_id_S) + "compL_{}".format(global_comp_idx)
            model_parameter_dict, _, _ = get_psbl_lightcurve_parameters(event_table_row, comp_table_rows, comp_idx, photometric_system, filter_name, event_id = 0)
            model = psbl_model_gen(model_parameter_dict)
            param_dict = lightcurve_parameter_gen(model, model_parameter_dict, np.array([global_comp_idx]), obj_id_L, obj_id_S, name, save_phot, phot_dir, overwrite)
            lightcurve_dict = {'obj_id_L' : obj_id_L, 'obj_id_S' : obj_id_S, 'companion_id_L' : comp_table_rows['companion_idx'][comp_idx], 
                                   'companion_id_S' : np.nan, 'class' : event_type}
            lightcurve_parameters.append([param_dict, lightcurve_dict])
            
            del global_comp_idx, param_dict,
            lightcurve_dict, model, model_parameter_dict, name
    
    elif event_type == 'BSPL':
        for comp_idx in range(len(comp_table_rows)):
            global_comp_idx = comp_table_rows['companion_idx'][comp_idx]
            name = "L_{}_S_{}".format(obj_id_L, obj_id_S) + "compS_{}".format(global_comp_idx)
            model_parameter_dict, _, _ = get_bspl_lightcurve_parameters(event_table_row, comp_table_rows, comp_idx, photometric_system, filter_name, red_law, event_id = 0)
            model = bspl_model_gen(model_parameter_dict)
            param_dict = lightcurve_parameter_gen(model, model_parameter_dict, np.array([global_comp_idx]), obj_id_L, obj_id_S, name, save_phot, phot_dir, overwrite)
            lightcurve_dict = {'obj_id_L' : obj_id_L, 'obj_id_S' : obj_id_S, 'companion_id_L' : np.nan, 
                                   'companion_id_S' : comp_table_rows['companion_idx'][comp_idx], 'class' : event_type}
            lightcurve_parameters.append([param_dict, lightcurve_dict])
            
            del global_comp_idx, param_dict,
            lightcurve_dict, model, model_parameter_dict, name

    # If multiple lightcurves were generated (i.e. for triples)
    # then take the lightcurve with the highest delta_m.
    # If all the delta_ms are nan, take the first one
    if len(lightcurve_parameters) == 1:
        lightcurve_parameters[0][0]['used_lightcurve'] = 1
    else:
        delta_ms = [x[0]['bin_delta_m'] for x in lightcurve_parameters]
        # if delta m is nan for all lightcurves, take the first one arbitrarily
        if np.all(np.isnan(delta_ms)):
            lightcurve_parameters[0][0]['used_lightcurve'] = 1
        else:
            #argmax ignoring nans
            max_delta_m_idx = np.nanargmax(delta_ms)
            lightcurve_parameters[max_delta_m_idx][0]['used_lightcurve'] = 1
            
            del max_delta_m_idx
        del delta_ms
    
    del multi_L, multi_S
    
    return lightcurve_parameters

        
def model_param_dict2fits_header(model_parameter_dict, phot_dir, name):
    """
    If photometry is being saved, saves the parameters
    of the model to generate the lightcurve in the first header
    
    Parameters
    ----------
    model_parameter_dict : dictionary
        Dictionary of the bagel model parameters 
        
    phot_dir : str
        Name of the directory photometry
    
    name : str
        Base name of fits file of photometry without the _phot 
        (i.e. L_{obj_id_L}_S_{obj_id_S}).
    """
    
    fits_file = phot_dir + '/' + name + '_phot.fits'
    with fits.open(fits_file, 'update') as f:
        hdr = f[0].header
        for key in list(model_parameter_dict.keys()):
            try:
                hdr[key] = model_parameter_dict[key]
            except ValueError:
                hdr[key] = str(model_parameter_dict[key])

    return

def lightcurve_parameter_gen(model, model_parameter_dict, comp_idxs, obj_id_L, obj_id_S,
                             name=None, save_phot=False, phot_dir=None, overwrite=False):
    """
    Find the parameters
    
    Parameters
    ----------
    model : function
        Function associated with generating the microlensing model (must be same as parameter_dict)
     
    model_parameter_dict : dict
        Dictionary of the bagel model parameters
    
    comp_idx : list
        List of indices into the comp_table of the companion(s) for which the
        binary event is being calculated.
    
    obj_id_L : int
        Object id of the lens associated with event
    
    obj_id_S : int
        Object id of the source associated with event
    
    name : str or None, optional
        Name of fits file to be saved.
        Default is None.
        
    save_phot : bool, optional
        If set to True, saves the photometry generated instead of just parameters.
        Default is False.
    
    phot_dir : str or None, optional
        Name of the directory photometry is saved if save_phot = True.
        This parameters is NOT optional if save_phot = True.
        Default is None.
        
    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.
    
    Returns
    -------
    param_dict : dict
        Dictionary of additional parameters added to companion table
        and one called 'mp_rows' which are the rows to be added to the multi peak
        table.
    """
    
    
    param_dict = {'n_peaks' : 0, 'bin_delta_m' : np.nan, 'tE_sys' : np.nan, 
                  'tE_primary' : np.nan, 'primary_t' : np.nan, 'avg_t' : np.nan, 
                  'std_t' : np.nan, 'asymmetry' : np.nan, 'mp_rows' : [], 'used_lightcurve' : 0}

    # Handles the case of modeling a triple source as a binary source,
    # But it's CO + CO + star and you're modeling the CO + CO pair (so no flux)
    if 'mag_src_sec' in list(model_parameter_dict.keys()):
        if np.isnan(model_parameter_dict['mag_src_sec']) and np.isnan(model_parameter_dict['mag_src_pri']):
            return param_dict
        
    # These covers extreme cases that will never
    # be detectable with such far t0s that jplephem
    # does not cover their date range
    jplephem_lower_mjd = 2287184.5 - 2400000.5
    jplephem_upper_mjd = 2688976.5 - 2400000.5
    if model.t0 < jplephem_lower_mjd or model.t0 > jplephem_upper_mjd:
        return param_dict
    
    # Calculate the photometry 
    duration=1000 # days
    time_steps=5000
    tmin = model.t0 - (duration / 2.0)
    tmax = model.t0 + (duration / 2.0)
    dt = np.linspace(tmin, tmax, time_steps)
    
    if 'get_all_arrays' in dir(model):
        img, amp = model.get_all_arrays(dt)
        phot = model.get_photometry(dt, amp_arr=amp)
    else: #BSPL
        phot = model.get_photometry(dt)

    if save_phot == True:
        if not os.path.exists(phot_dir):
            os.makedirs(phot_dir)
        foo = Table((dt, phot), names=['time', 'phot'])
        foo.write(phot_dir + '/' + name + '_phot.fits', overwrite=overwrite)
        model_param_dict2fits_header(model_parameter_dict, phot_dir, name)

    #because this is magnitudes max(phot) is baseline and min(phot) is peak
    #baseline 2000tE away get_photometry
    bin_delta_m = max(phot) - min(phot)
    param_dict['bin_delta_m'] = bin_delta_m

    # system tE is the time between where the event is first above 10% of it's base mag
    # to the last time. This may not happen for some events
    tenp = np.where(phot < (max(phot) - 0.1*bin_delta_m))[0]
    if len(tenp) == 0:
        param_dict['tE_sys'] = np.nan
    else:
        param_dict['tE_sys'] = max(dt[tenp]) - min(dt[tenp])

    # Find peaks
    peaks, _ = find_peaks(-phot, prominence = 10e-5, width =1) 
    param_dict['n_peaks'] = len(peaks)
    if len(peaks) == 0:
        return param_dict
        
    param_dict['primary_t'] = dt[peaks][np.argmin(phot[peaks])]
    param_dict['avg_t'] = np.average(dt[peaks])
    param_dict['std_t'] = np.std(dt[peaks])
    
    # Find asymmetry (0 if symmetric, larger if asymmetric)
    # Uses 50 degree chebyshev polynomial (see eq 7 in Night et al. 2010)
    indices = np.arange(0,51)
    cheb_fit = np.polynomial.chebyshev.Chebyshev.fit(dt, phot, 50).coef
    odd_cheb = cheb_fit[indices%2==1]
    even_cheb = cheb_fit[indices%2==0]
    asymm = np.sqrt(np.sum(odd_cheb**2)/np.sum(even_cheb**2))
    param_dict['asymmetry'] = asymm
    
    # Split up peaks by minima between peaks
    # Note since it's magnitudes all the np.min and such are maxima
    split_data = []
    start_idx = 0
    if len(peaks) > 1:
        for i in range(len(peaks) - 1):
            min_btwn_peaks = np.max(phot[peaks[i]:peaks[i+1]])
            end_idx = np.where(phot[start_idx:] == min_btwn_peaks)[0][0] + start_idx
            split_data.append([dt[start_idx:end_idx], phot[start_idx:end_idx]])
            start_idx = end_idx
    split_data.append([dt[start_idx:], phot[start_idx:]])
    split_data = np.array(split_data, dtype='object')
    
    highest_peak = np.argmin(phot[peaks])

    if len(split_data[highest_peak][1]) != 0:
        highest_bump_mag = max(split_data[highest_peak][1]) - min(split_data[highest_peak][1])
        highest_half = np.where(split_data[highest_peak][1] < (max(split_data[highest_peak][1]) - 0.5*highest_bump_mag))[0]
        
        if len(highest_half) != 0:
            tE_primary = max(dt[highest_half]) - min(dt[highest_half])
            param_dict['tE_primary'] = tE_primary
            
    # For events with more than one peak, add them to the multi peak table
    rows = []
    if len(peaks) > 1:
        n_peaks = len(peaks)
        for i in range(len(peaks)):
            t = dt[peaks[i]]
            delta_m = max(phot) - phot[peaks[i]]
            ratio = np.min(phot[peaks])/phot[peaks[i]]

            # Don't log primary peak
            if np.min(phot[peaks]) == phot[peaks[i]]:
                continue

            if len(split_data[i][1]) == 0:
                tE = np.nan
                rows.append([comp_idxs, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
                continue

            split_bump_mag = max(split_data[i][1]) - min(split_data[i][1])
            split_half = np.where(split_data[i][1] < (max(split_data[i][1]) - 0.5*split_bump_mag))[0]
            if len(split_half) == 0:
                tE = np.nan
                rows.append([comp_idxs, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
                continue

            tE = max(dt[split_half]) - min(dt[split_half])
            rows.append([comp_idxs, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
            
        param_dict['mp_rows'] = rows
        
    return param_dict

def get_psbl_lightcurve_parameters(event_table, comp_table, comp_idx, photometric_system, filter_name, event_id = None):
    """
    Find the parameters for PSBL_PhotAstrom_Par_Param7 from 
    event_table and comp_table.

    Parameters
    ----------
    event_table : Astropy table
        Table containing the events calculated from refine_events.
    
    comp_table : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx : int
        Index into the comp_table of the companion for which the psbl is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id.
        Default is None.
        
    Returns
    -------
    psbl_parameter_dict : dict
        Dictionary of the PSBL_PhotAstrom_Par_Param7 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event
        
    """
    obj_id_L = comp_table['obj_id_L'][comp_idx]
    obj_id_S = comp_table['obj_id_S'][comp_idx]
    
    if event_id is None:
        event_id = (np.where(np.logical_and((event_table['obj_id_L'] == obj_id_L), (event_table['obj_id_S'] == obj_id_S)))[0])[0]
        
    L_coords = SkyCoord(l = event_table[event_id]['glon_L']*unit.degree, b = event_table[event_id]['glat_L']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_L']*unit.mas/unit.year, 
                            pm_b = event_table[event_id]['mu_b_L']*unit.mas/unit.year, frame ='galactic')
    S_coords = SkyCoord(l = event_table[event_id]['glon_S']*unit.degree, b = event_table[event_id]['glat_S']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_S']*unit.mas/unit.year, 
                              pm_b = event_table[event_id]['mu_b_S']*unit.mas/unit.year, frame ='galactic')
    raL = L_coords.icrs.ra.value # Lens R.A.
    decL = L_coords.icrs.dec.value # Lens dec
    mL1 = event_table[event_id]['mass_L'] # msun (Primary lens current mass)
    mL2 = comp_table['mass'][comp_idx] # msun (Companion lens current mass)
    t0 = event_table[event_id]['t0'] # mjd
    xS0 = np.array([0, 0]) #arbitrary offset (arcsec)
    beta = event_table[event_id]['u0']*event_table[event_id]['theta_E']#5.0
    muL = np.array([L_coords.icrs.pm_ra_cosdec.value, L_coords.icrs.pm_dec.value]) #lens proper motion mas/year
    muS = np.array([S_coords.icrs.pm_ra_cosdec.value, S_coords.icrs.pm_dec.value]) #source proper motion mas/year
    dL = event_table[event_id]['rad_L']*10**3 #Distance to lens
    dS = event_table[event_id]['rad_S']*10**3 #Distance to source
    sep = comp_table['sep'][comp_idx] #mas (separation between primary and companion)
    alpha = comp_table['alpha'][comp_idx]
    mag_src = event_table[event_id]['%s_%s_app_S' % (photometric_system, filter_name)]
    b_sff = event_table[event_id]['f_blend_%s' % filter_name] #ASSUMES ALL BINARY LENSES ARE BLENDED
    model_name = 'PSBL_PhotAstrom_Par_Param7'
    
    psbl_parameter_dict = {'raL': raL, 'decL': decL, 'mL1': mL1, 'mL2': mL2, 
                           't0': t0, 'xS0': xS0, 'beta': beta, 'muL': muL, 
                           'muS': muS, 'dL': dL, 'dS': dS, 'sep': sep, 
                           'alpha': alpha, 'mag_src': mag_src, 'b_sff': b_sff, 
                           'model': model_name}
    return psbl_parameter_dict, obj_id_L, obj_id_S
    

def psbl_model_gen(psbl_parameter_dict):
    """
    Generate psbl_photastrom_par_param1 model from parameter dict
    
    Parameters
    ----------
    psbl_parameter_dict : dict
        Dictionary of the PSBL_PhotAstrom_Par_Param7 parameters 
        
    Returns
    -------
    psbl model
    """
    ##########
    # Calculate binary model and photometry
    ##########
    raL = psbl_parameter_dict['raL'] # Lens R.A.
    decL = psbl_parameter_dict['decL'] # Lens dec
    mL1 = psbl_parameter_dict['mL1'] # msun (Primary lens current mass)
    mL2 = psbl_parameter_dict['mL2'] # msun (Companion lens current mass)
    t0 = psbl_parameter_dict['t0'] # mjd
    xS0 = psbl_parameter_dict['xS0'] #arbitrary offset (arcsec)
    beta = psbl_parameter_dict['beta']
    muL = psbl_parameter_dict['muL'] #lens proper motion mas/year
    muS = psbl_parameter_dict['muS'] #source proper motion mas/year
    dL = psbl_parameter_dict['dL'] #Distance to lens
    dS = psbl_parameter_dict['dS'] #Distance to source
    sep = psbl_parameter_dict['sep'] #mas (separation between primary and companion)
    alpha = psbl_parameter_dict['alpha']
    mag_src = psbl_parameter_dict['mag_src']
    b_sff = psbl_parameter_dict['b_sff'] #ASSUMES ALL BINARY LENSES ARE BLENDED

    psbl = model.PSBL_PhotAstrom_Par_Param7(mL1, mL2, t0, xS0[0], xS0[1],
                               beta, muL[0], muL[1], muS[0], muS[1], dL, dS,
                               sep, alpha, [b_sff], [mag_src], 
                               raL=raL, decL=decL, 
                               root_tol = 0.00000001)
    return psbl


def get_bspl_lightcurve_parameters(event_table, comp_table, comp_idx, photometric_system, filter_name, red_law, event_id = None):
    """
    Find the parameters for BSPL_PhotAstrom_Par_Param1 from 
    event_table and comp_table.

    Parameters
    ----------
    event_table : Astropy table
        Table containing the events calculated from refine_events.
    
    comp_table : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx : int
        Index into the comp_table of the companion for which the psbl is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    red_law : str
        Redenning law
    
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id
        
    Returns
    -------
    bspl_parameter_dict : dict
        Dictionary of the BSPL_PhotAstrom_Par_Param1 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event
        
    """
    obj_id_L = comp_table['obj_id_L'][comp_idx]
    obj_id_S = comp_table['obj_id_S'][comp_idx]
    
    if event_id is None:
        event_id = (np.where(np.logical_and((event_table['obj_id_L'] == obj_id_L), (event_table['obj_id_S'] == obj_id_S)))[0])[0]
    L_coords = SkyCoord(l = event_table[event_id]['glon_L']*unit.degree, b = event_table[event_id]['glat_L']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_L']*unit.mas/unit.year, 
                            pm_b = event_table[event_id]['mu_b_L']*unit.mas/unit.year, frame ='galactic')
    S_coords = SkyCoord(l = event_table[event_id]['glon_S']*unit.degree, b = event_table[event_id]['glat_S']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_S']*unit.mas/unit.year, 
                              pm_b = event_table[event_id]['mu_b_S']*unit.mas/unit.year, frame ='galactic')
    f_i = filt_dict[photometric_system + '_' + filter_name][red_law]
    abs_mag_sec = comp_table['m_%s_%s' % (photometric_system, filter_name)][comp_idx]
    
    raL = L_coords.icrs.ra.value # Lens R.A.
    decL = L_coords.icrs.dec.value # Lens dec
    mL = event_table[event_id]['mass_L'] # msun (Lens current mass)
    t0 = event_table[event_id]['t0'] # mjd
    beta = event_table[event_id]['u0']*event_table[event_id]['theta_E']#5.0
    dL = event_table[event_id]['rad_L']*10**3 #Distance to lens
    dL_dS = dL/(event_table[event_id]['rad_S']*10**3) #Distance to lens/Distance to source
    xS0 = np.array([0, 0]) #arbitrary offset (arcsec)
    muL_E = L_coords.icrs.pm_ra_cosdec.value #lens proper motion mas/year
    muL_N = L_coords.icrs.pm_dec.value #lens proper motion mas/year
    muS_E = S_coords.icrs.pm_ra_cosdec.value #lens proper motion mas/year
    muS_N = S_coords.icrs.pm_dec.value #lens proper motion mas/year
    sep = comp_table['sep'][comp_idx] #mas (separation between primary and companion)
    alpha = comp_table['alpha'][comp_idx]
    mag_src_sec = calc_app_mag(event_table[event_id]['rad_S'], abs_mag_sec, event_table[event_id]['exbv_S'], f_i)
    mag_src_pri = binary_utils.subtract_magnitudes(event_table[event_id]['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec)
    b_sff = event_table[event_id]['f_blend_%s' % filter_name] #ASSUMES THAT SOURCE BINARIES ARE BLENDED
    model_name = 'BSPL_PhotAstrom_Par_Param1'
    
    bspl_parameter_dict = {'model': model_name, 'raL': raL, 'decL': decL, 'mL': mL,
                           't0': t0, 'xS0': xS0, 'beta': beta, 
                           'muL_E': muL_E, 'muL_N': muL_N, 'muS_E': muS_E, 'muS_N': muS_N,
                           'dL': dL, 'dL_dS': dL_dS, 'sep': sep, 
                           'alpha': alpha, 'mag_src_pri': mag_src_pri, 'mag_src_sec': mag_src_sec, 
                           'b_sff': b_sff}
    
    return bspl_parameter_dict, obj_id_L, obj_id_S

def bspl_model_gen(bspl_parameter_dict):
    """
    Generate bspl_photastrom_par_param1 model from parameter dict
    
    Parameters
    ----------
    bspl_parameter_dict : dict
        Dictionary of the BSPL_PhotAstrom_Par_Param1 parameters 
        
    Returns
    -------
    bspl model
    """
    ##########
    # Calculate binary model and photometry
    ##########
    raL = bspl_parameter_dict['raL'] # Lens R.A.
    decL = bspl_parameter_dict['decL'] # Lens dec
    mL = bspl_parameter_dict['mL'] # msun (Lens current mass)
    t0 = bspl_parameter_dict['t0'] # mjd
    beta = bspl_parameter_dict['beta']
    dL = bspl_parameter_dict['dL'] #Distance to lens
    dL_dS = bspl_parameter_dict['dL_dS'] #Distance to lens/Distance to source
    xS0 = bspl_parameter_dict['xS0'] #arbitrary offset (arcsec)
    muL_E = bspl_parameter_dict['muL_E'] #lens proper motion mas/year
    muL_N = bspl_parameter_dict['muL_N'] #lens proper motion mas/year
    muS_E = bspl_parameter_dict['muS_E'] #lens proper motion mas/year
    muS_N = bspl_parameter_dict['muS_N'] #lens proper motion mas/year
    sep = bspl_parameter_dict['sep'] #mas (separation between primary and companion)
    alpha = bspl_parameter_dict['alpha']
    mag_src_sec = bspl_parameter_dict['mag_src_sec']
    mag_src_pri = bspl_parameter_dict['mag_src_pri']
    b_sff = bspl_parameter_dict['b_sff'] #ASSUMES ALL BINARY LENSES ARE BLENDED

    bspl = model.BSPL_PhotAstrom_Par_Param1(mL, t0, beta, dL, dL_dS, 
                                   xS0[0], xS0[1], muL_E, muL_N, muS_E, muS_N,
                                   sep, alpha, [mag_src_pri], [mag_src_sec], [b_sff],
                                   raL=raL, decL=decL)
    
    return bspl

def get_bsbl_lightcurve_parameters(event_table, comp_table, comp_idx_L, comp_idx_S, photometric_system, filter_name, red_law, event_id = None):
    """
    Find the parameters for BSPL_PhotAstrom_Par_Param2 from 
    event_table and comp_table.

    Parameters
    ----------
    event_table : Astropy table
        Table containing the events calculated from refine_events.
    
    comp_table : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx_L : int
        Index into the comp_table of the lens companion for which the model is being calculated.
        
    comp_idx_S : int
        Index into the comp_table of the source companion for which the model is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    red_law : str
        Redenning law
        
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id.
        Default is None.
        
    Returns
    -------
    bsbl_parameter_dict : dict
        Dictionary of the BSBL_PhotAstrom_Par_Param2 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event
        
    """
    obj_id_L = comp_table['obj_id_L'][comp_idx_L] # This is equivalent to doing comp_idx_S
    obj_id_S = comp_table['obj_id_S'][comp_idx_S]

    if event_id is None:
        event_id = (np.where(np.logical_and((event_table['obj_id_L'] == obj_id_L), (event_table['obj_id_S'] == obj_id_S)))[0])[0]
    
    L_coords = SkyCoord(l = event_table[event_id]['glon_L']*unit.degree, b = event_table[event_id]['glat_L']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_L']*unit.mas/unit.year, 
                            pm_b = event_table[event_id]['mu_b_L']*unit.mas/unit.year, frame ='galactic')    
    S_coords = SkyCoord(l = event_table[event_id]['glon_S']*unit.degree, b = event_table[event_id]['glat_S']*unit.degree, 
                            pm_l_cosb = event_table[event_id]['mu_lcosb_S']*unit.mas/unit.year, 
                              pm_b = event_table[event_id]['mu_b_S']*unit.mas/unit.year, frame ='galactic')

    f_i = filt_dict[photometric_system + '_' + filter_name][red_law]
    abs_mag_sec = comp_table['m_%s_%s' % (photometric_system, filter_name)][comp_idx_S]

    raL = L_coords.icrs.ra.value # Lens R.A.
    decL = L_coords.icrs.dec.value # Lens dec
    mLp = event_table[event_id]['mass_L'] # msun (Lens current mass)
    mLs = comp_table['mass'][comp_idx_L] # msun (Companion lens current mass)
    t0 = event_table[event_id]['t0'] # mjd
    beta = event_table[event_id]['u0']*event_table[event_id]['theta_E']#5.0
    dL = event_table[event_id]['rad_L']*10**3 #Distance to lens
    dS = event_table[event_id]['rad_S']*10**3 #Distance to source
    xS0_E = 0.0 #arbitrary offset (arcsec)
    xS0_N = 0.0 #arbitrary offset (arcsec)
    muL_E = L_coords.icrs.pm_ra_cosdec.value #lens proper motion mas/year
    muL_N = L_coords.icrs.pm_dec.value #lens proper motion mas/year
    muS_E = S_coords.icrs.pm_ra_cosdec.value #lens proper motion mas/year
    muS_N = S_coords.icrs.pm_dec.value #lens proper motion mas/year
    sepL = comp_table['sep'][comp_idx_L] #mas (separation between primary and companion)
    alphaL = comp_table['alpha'][comp_idx_L] # PA of binary on the sky
    sepS = comp_table['sep'][comp_idx_S] #mas (separation between primary and companion)
    alphaS = comp_table['alpha'][comp_idx_S] # PA of source binary on the sky
    mag_src_sec = calc_app_mag(event_table[event_id]['rad_S'], abs_mag_sec, event_table[event_id]['exbv_S'], f_i)
    mag_src_pri = binary_utils.subtract_magnitudes(event_table[event_id]['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec)
    b_sff = event_table[event_id]['f_blend_%s' % filter_name] #ASSUMES THAT SOURCE BINARIES ARE BLENDED
    model_name = 'BSBL_PhotAstrom_Par_Param2'
    
    bsbl_parameter_dict = {'model': model_name, 'raL': raL, 'decL': decL, 'mLp': mLp, 'mLs': mLs,
                           't0': t0, 'xS0_E': xS0_E, 'xS0_N': xS0_N, 'beta': beta, 
                           'muL_E': muL_E, 'muL_N': muL_N, 'muS_E': muS_E, 'muS_N': muS_N,
                           'dL': dL, 'dS': dS, 'sepL': sepL, 'alphaL': alphaL, 
                           'sepS': sepS, 'alphaS': alphaS,
                           'mag_src_pri': mag_src_pri, 'mag_src_sec': mag_src_sec, 
                           'b_sff': b_sff}
        
    return bsbl_parameter_dict, obj_id_L, obj_id_S

def bsbl_model_gen(bsbl_parameter_dict):
    """
    Generate bsbl_photastrom_par_param1 model from parameter dict
    
    Parameters
    ----------
    bspl_parameter_dict : dict
        Dictionary of the BSBL_PhotAstrom_Par_Param1 parameters 
        
    Returns
    -------
    bsbl model
    """
    ##########
    # Calculate binary model and photometry
    ##########
    raL = bsbl_parameter_dict['raL'] # Lens R.A.
    decL = bsbl_parameter_dict['decL'] # Lens dec
    mLp = bsbl_parameter_dict['mLp'] # msun (Lens current mass)
    mLs = bsbl_parameter_dict['mLs'] # msun (Lens companion current mass)
    t0 = bsbl_parameter_dict['t0'] # mjd
    beta = bsbl_parameter_dict['beta']
    dL = bsbl_parameter_dict['dL'] #Distance to lens
    dS = bsbl_parameter_dict['dS'] #Distance to source
    xS0_E = bsbl_parameter_dict['xS0_E'] #arbitrary offset (arcsec)
    xS0_N = bsbl_parameter_dict['xS0_N'] #arbitrary offset (arcsec)
    muL_E = bsbl_parameter_dict['muL_E'] #lens proper motion mas/year
    muL_N = bsbl_parameter_dict['muL_N'] #lens proper motion mas/year
    muS_E = bsbl_parameter_dict['muS_E'] #lens proper motion mas/year
    muS_N = bsbl_parameter_dict['muS_N'] #lens proper motion mas/year
    sepL = bsbl_parameter_dict['sepL'] #mas (separation between primary and companion)
    alphaL = bsbl_parameter_dict['alphaL']
    sepS = bsbl_parameter_dict['sepS'] #mas (separation between primary and companion)
    alphaS = bsbl_parameter_dict['alphaS']
    mag_src_sec = bsbl_parameter_dict['mag_src_sec']
    mag_src_pri = bsbl_parameter_dict['mag_src_pri']
    b_sff = bsbl_parameter_dict['b_sff'] #ASSUMES ALL BINARY LENSES ARE BLENDED

    bsbl = model.BSBL_PhotAstrom_Par_Param2(mLp, mLs, t0, xS0_E, xS0_N,
                                              beta, muL_E, muL_N, muS_E, muS_N,
                                              dL, dS, sepL, alphaL, sepS, alphaS,
                                              mag_src_pri, mag_src_sec, b_sff,
                                              raL=raL, decL=decL,
                                              root_tol=1e-4)
    return bsbl

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
        if ii <= 1 or ii >= len(log_list) - 2:
            continue

        log_list_entry = log_list[ii].replace(b'# ',b'').split()

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
        <h5file_name>_label.fits (after removing the \*.h5).

    Returns
    -------
    label_file : Astropy table
        Table containing the dataset name, and corresponding l, b, and number of objects.

    """
    label_file_start_time = time.time()
    
    hf = h5py.File(h5file_name + '.h5', 'r')
    l_array = hf['long_bin_edges']
    b_array = hf['lat_bin_edges']
    
    dset_num = (len(l_array)-1)*(len(b_array)-1)
    
    data_dict = {'file_name': np.chararray(dset_num,itemsize=20), 'long_start': np.empty(dset_num), 'long_end': np.empty(dset_num),
                 'lat_start': np.empty(dset_num), 'lat_end': np.empty(dset_num), 'objects': np.empty(dset_num, dtype = int),
                 'N_stars': np.empty(dset_num, dtype = int), 'N_WD': np.empty(dset_num, dtype = int),
                 'N_NS': np.empty(dset_num, dtype = int), 'N_BH': np.empty(dset_num, dtype = int)}
    
    dset_num = 0
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
                values, counts = np.unique(dataset['rem_id'], return_counts = True)
                if 0 in values:
                    N_stars = counts[values == 0][0]
                else:
                    N_stars = 0
                if 101 in values:
                    N_WD = counts[values == 101][0]
                else:
                    N_WD = 0
                if 102 in values:
                    N_NS = counts[values == 102][0]
                else:
                    N_NS = 0
                if 103 in values:
                    N_BH = counts[values == 103][0]
                else:
                    N_BH = 0
            
            data_dict['file_name'][dset_num] = dset_name
            data_dict['long_start'][dset_num] = l_array[ll]
            data_dict['long_end'][dset_num] = l_array[ll + 1]
            data_dict['lat_start'][dset_num] = b_array[bb]
            data_dict['lat_end'][dset_num] = b_array[bb + 1]
            data_dict['objects'][dset_num] = len(dataset)
            data_dict['N_stars'][dset_num] = N_stars
            data_dict['N_WD'][dset_num] = N_WD
            data_dict['N_NS'][dset_num] = N_NS
            data_dict['N_BH'][dset_num] = N_BH
            
            dset_num += 1
            del dataset

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
    
    print('make label file time: {0:f} s'.format(time.time() - label_file_start_time))

    return

def make_label_file_no_bins(h5file_name, overwrite=False):
    """
    Given an output root for an .h5 file, creates a table of
    dataset name and number of objects.
    Writes out the Astropy table as a .fits file.
    Used when binning = False (when simulating full sky)

    Parameters
    ----------
    h5file_name : string
        The path and name of the input h5file containing the
        population synthesis results. We will read this in and
        make a new output file entitled:
        <h5file_name>_label.fits (after removing the \*.h5).

    Returns
    -------
    label_file : Astropy table
        Table containing the dataset name, and corresponding l, b, and number of objects.

    """
    data_dict = {'file_name': [], 'objects': [],
                 'N_stars': [], 'N_WD': [], 'N_NS': [], 'N_BH': []}

    hf = h5py.File(h5file_name + '.h5', 'r')

    dataset = hf['objects']

    if len(dataset) == 0:
        N_stars = 0
        N_WD = 0
        N_NS = 0
        N_BH = 0
    else:
        values, counts = np.unique(dataset['rem_id'], return_counts = True)
        if 0 in values:
            N_stars = counts[values == 0][0]
        else:
            N_stars = 0
        if 101 in values:
            N_WD = counts[values == 101][0]
        else:
            N_WD = 0
        if 102 in values:
            N_NS = counts[values == 102][0]
        else:
            N_NS = 0
        if 103 in values:
            N_BH = counts[values == 103][0]
        else:
            N_BH = 0

    data_dict['file_name'].append('object')
    data_dict['objects'].append(dataset.shape[0])
    data_dict['N_stars'].append(N_stars)
    data_dict['N_WD'].append(N_WD)
    data_dict['N_NS'].append(N_NS)
    data_dict['N_BH'].append(N_BH)

    hf.close()

    label_file = Table(data_dict, names=('file_name', 'objects',
                                         'N_stars', 'N_WD', 'N_NS', 'N_BH'))

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
    ----------
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
        Lens current mass, in solar masses

    d_L : float
        Distance to lens, in kpc

    d_S : float
        Distance to source, in kpc

    Returns
    -------
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

    Returns
    -------
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




def calc_magnification(u):
    """
    Calculate the magnification factor A(u)

    Parameters
    ----------
    u : float or array
        Dimensionless lens-source separation, normalized in Einstein radius units

    Returns
    -------
    Magnification factor : float or array

    """
    u = np.abs(u)
    return (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))


def calc_delta_c(u, thetaE):
    """
    Calculate the maximum centroid shift for a dark lens,
    no neighbors

    Parameters
    ----------
    u : float or array
        Dimensionless lens-source separation, normalized in Einstein radius units.

    thetaE : float
        Einstein radius of lens.

    Returns
    -------
    delta_c : float or array
        Maximum centroid shift
    """
    u = np.abs(u)
    delta_c = u * thetaE / (u ** 2 + 2)

    return delta_c


def calc_bump_amp(u0, f_S, f_L, f_N):
    """
    Calculate the "bump" amplitude, given the minimum separation and the fluxes
    of the (unmagnified) source, lens, and neighbors.
    The bump amplitude (mags) is:
    abs(m_peak - m_base) = -2.5 log10 ((A(u0) * f_S + f_L + f_N)/(f_S + f_L + f_N))

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

    Returns
    -------
    m_bump : float or array
        Bump magnitude
    """
    A = calc_magnification(np.abs(u0))
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

    Returns
    -------
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

    Returns
    -------
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

    Returns
    -------
    m_E : float or array
        Magnitude of extinction

    """
    m_E = E * f

    return m_E


def get_Alambda_AKs(red_law_name, lambda_eff):
    """
    Get Alambda/AKs. NOTE: this doesn't work for every law in SPISEA!
    Naming convention is not consistent. Change SPISEA or add if statements?

    Parameters
    ----------
    red_law_name : str
        The name of the reddening law
    lambda_eff : float
        Wavelength in microns

    Returns
    -------
    Alambda_AKs : float
        Alambda/AKs

    """
    red_law_class = getattr(reddening, 'RedLaw' + red_law_name)
    red_law = red_law_class()
    red_law_method = getattr(red_law, red_law_name)
    Alambda_AKs = red_law_method(lambda_eff, 1)

    return Alambda_AKs

