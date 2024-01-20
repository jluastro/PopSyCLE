#! /usr/bin/env python
"""
run.py
Executable to run the PopSyCLE pipeline.
"""
import inspect
import os
from pathlib import Path
import argparse
from argparse import RawTextHelpFormatter
import yaml
import sys
import time
import glob
from popsycle import synthetic
from popsycle import utils
from popsycle.synthetic import _check_run_galaxia
from popsycle.synthetic import _check_perform_pop_syn
from popsycle.synthetic import _check_calc_events
from popsycle.synthetic import _check_refine_events
from popsycle.synthetic import _check_refine_binary_events
from popsycle.synthetic import multiplicity_list


def _return_filename_dict(output_root):
    """
    Return the filenames of the files output by the pipeline

    Parameters
    ----------
    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    Output
    ------
    filename_dict : dict
        Dictionary containing the names of the files output by the pipeline

    """
    # Write out all of the filenames using the output_root
    ebf_filename = '%s.ebf' % output_root
    hdf5_filename = '%s.h5' % output_root
    events_filename = '%s_events.fits' % output_root
    blends_filename = '%s_blends.fits' % output_root
    noevents_filename = '%s_NOEVENTS.txt' % output_root

    # Add the filenames to a dictionary
    filename_dict = {
        'ebf_filename': ebf_filename,
        'hdf5_filename': hdf5_filename,
        'events_filename': events_filename,
        'blends_filename': blends_filename,
        'noevents_filename': noevents_filename
    }

    return filename_dict


def _check_for_output(filename, overwrite=False):
    """
    Checks for the existence of files and either overwrites them or
    raises a warning.

    Parameters
    ----------
    filename : str
        Name of the file to be inspected

    overwrite : bool
        Flag to determine whether to overwrite in the presence of the file
        or to raise an error. If True, file is overwritten.
        If False, error is raised. Default False.

    Output
    ------
    status : bool
        Status of operation.
        True: Error due to already existing file.
        False: File either does not exist or was successfully deleted.

    """
    if os.path.exists(filename):
        if overwrite:
            os.remove(filename)
            return False
        else:
            print('Error: Output {0} exists, cannot continue. Either '
                  'rename {0} or rerun run.py with the --overwrite '
                  'flag.'.format(filename))
            return True
    else:
        return False


def generate_field_config_file(longitude, latitude, area,
                               config_filename='field_config.yaml'):
    """
    Save field configuration parameters from a dictionary into a yaml file

    Parameters
    ----------
    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    Optional Parameters
    -------------------
    config_filename : str
        Name of the configuration file
        Default: field_config.yaml

    Output
    ------
    None
    """

    config = {'longitude': longitude,
              'latitude': latitude,
              'area': area}
    generate_config_file(config_filename, config)


def generate_slurm_config_file(path_python='python', account='ulens',
                               queue='regular', resource='cpu',
                               memory=512, n_cores_per_node=64, n_nodes_max=3072,
                               memory_max=512,
                               walltime_max='12:00:00',
                               additional_lines=['module load cray-hdf5/1.10.5.2',
                                                 'export HDF5_USE_FILE_LOCKING=FALSE'],
                               config_filename='slurm_config.yaml'):
    """
    Save slurm configuration parameters from a dictionary into a yaml file

    Parameters
    ----------
    path_python : str
        Path to the python executable

    account : str
        Project account name to charge

    queue : str
        Scheduler queue type

    resource : str
        Computing resource name

    memory : int
        Amount of memory allocated for the job, in units of GB

    n_cores_per_node : int
        Number of cores in each node of the compute resource

    n_nodes_max : int
        Total number of nodes in the compute resource

    memory_max : int
        Memory per node in the computer resource, in units of GB

    walltime_max : int
        Maximum number of hours for single job on the compute resource
        Format: hh:mm:ss

    additional_lines : list of strings
        Additional lines to be run before executing run.py

    Optional Parameters
    -------------------
    config_filename : str
        Name of the configuration file
        Default: slurm_config.yaml

    Output
    ------
    None
    """

    config = {'path_python': path_python,
              'account': account,
              'queue': queue,
              'resource': resource,
              'memory': memory,
              'additional_lines': additional_lines,
              resource: {'n_cores_per_node': n_cores_per_node,
                         'n_nodes_max': n_nodes_max,
                         'memory_max': memory_max,
                         'walltime_max': walltime_max}}
    generate_config_file(config_filename, config)


def generate_popsycle_config_file(radius_cut=2, obs_time=1000,
                                  n_obs=101, theta_frac=2, blend_rad=0.75,
                                  isochrones_dir='/Users/myself/popsycle_isochrones',
                                  IFMR='Raithel18',
                                  galaxia_galaxy_model_filename='/Users/myself/galaxia_galaxy_model_filename',
                                  bin_edges_number=None,
                                  BH_kick_speed_mean=50,
                                  NS_kick_speed_mean=400,
                                  photometric_system='ubv',
                                  filter_name='R', red_law='Damineli16',
                                  multiplicity=None,
                                  binning = True,
                                  config_filename='popsycle_config.yaml'):
    """
    Save popsycle configuration parameters from a dictionary into a yaml file

    Parameters
    ----------
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

    isochrones_dir : str
        Directory for SPISEA isochrones

    IFMR : string
        The name of the IFMR object from SPISEA. For more information on these objects see ifmr.py
        in SPISEA.
        'Raithel18' = IFMR_Raithel18
        'Spera15' = IFMR_Spera15
        'SukhboldN20' = IFMR_N20_Sukhbold

    galaxia_galaxy_model_filename : str
        Name of the galaxia galaxy model, as outlined at https://github.com/jluastro/galaxia

    bin_edges_number : int
        Number of edges for the bins
            bins = bin_edges_number - 1
        Total number of bins is
            N_bins = (bin_edges_number - 1)**2
        If None (default), then number of bins is
            bin_edges_number = int(60 * 2 * radius) + 1

    BH_kick_speed_mean : float
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.
        Defaults to 50 km/s.

    NS_kick_speed_mean : float
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.
        Defaults to 400 km/s based on distributions found by
        Hobbs et al 2005 'A statistical study of 233 pulsar proper motions'.
        https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract

    photometric_system : str
        The name of the photometric system in which the filter exists.

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    red_law : str
        The name of the reddening law to use from SPISEA.

    multiplicity: str
        If a resovled multiplicity object is specified,
        the table will be generated with resolved multiples.
        Default is None.
    
    binning : bool
        If set to True, bins files as specified by bin_edges_numbers or default.
        If set to False, no bins (SET TO FALSE IF DOING FULL SKY DOWNSAMPLED).
        Default is True.

    Optional Parameters
    -------------------
    config_filename : str
        Name of the configuration file
        Default: popsycle_config.yaml

    Output
    ------
    None
    """

    if bin_edges_number is None:
        bin_edges_number = 'None'
    if isochrones_dir == '/Users/myself/popsycle_isochrones':
        raise Exception("'isochrones_dir' must be set by the user. "
                        "The default value is only an example.")
    if galaxia_galaxy_model_filename == '/Users/myself/galaxia_galaxy_model_filename':
        raise Exception("'galaxia_galaxy_model_filename' must be set by the user. "
                        "The default value is only an example.")

    if multiplicity is None:
        multiplicity = 'None'
    if multiplicity not in multiplicity_list:
        raise Exception('multiplicity must be None or "ResolvedDK"')

    config = {'radius_cut': radius_cut,
              'obs_time': obs_time,
              'n_obs': n_obs,
              'theta_frac': theta_frac,
              'blend_rad': blend_rad,
              'isochrones_dir': os.path.abspath(isochrones_dir),
              'IFMR' : IFMR,
              'galaxia_galaxy_model_filename': os.path.abspath(galaxia_galaxy_model_filename),
              'bin_edges_number': bin_edges_number,
              'BH_kick_speed_mean': BH_kick_speed_mean,
              'NS_kick_speed_mean': NS_kick_speed_mean,
              'photometric_system': photometric_system,
              'filter_name': filter_name,
              'red_law': red_law,
              'multiplicity': multiplicity,
              'binning':binning}
    generate_config_file(config_filename, config)


def generate_config_file(config_filename, config):
    """
    Save configuration parameters from a dictionary into a yaml file

    Parameters
    ----------
    config_filename : str
        Name of the configuration file

    config : dict
        Dictionary containing the configuration parameters

    Output
    ------
    None

    """
    with open(config_filename, 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=True)


def load_config_file(config_filename):
    """
    Load configuration parameters from a yaml file into a dictionary

    Parameters
    ----------
    config_filename : str
        Name of the configuration file

    Output
    ------
    config : dict
        Dictionary containing the configuration parameters

    """
    with open(config_filename, 'r') as f:
        config = yaml.load(f, Loader=yaml.Loader)
    return config


def _check_slurm_config(slurm_config, walltime):
    """
    Checks that the values in slurm_config are valid

    Parameters
    ----------
    slurm_config : dict
        Dictionary of values for configuring slurm scripts

    walltime : str
        Amount of walltime that the script will request from slurm.
        Format: hh:mm:ss
    """

    if 'path_python' not in slurm_config:
        raise Exception('path_python must be set in slurm_config')

    path_python = slurm_config['path_python']
    if type('path_python') != str:
        raise Exception('path_python (%s) must be a string.' % str(path_python))

    if 'account' not in slurm_config:
        raise Exception('account must be set in slurm_config')

    account = slurm_config['account']
    if type('account') != str:
        raise Exception('account (%s) must be a string.' % str(account))

    if 'queue' not in slurm_config:
        raise Exception('queue must be set in slurm_config')

    queue = slurm_config['queue']
    if type('queue') != str:
        raise Exception('queue (%s) must be a string.' % str(queue))

    if 'resource' not in slurm_config:
        raise Exception('resource must be set in slurm_config')

    resource = slurm_config['resource']
    if type('resource') != str:
        raise Exception('resource (%s) must be a string.' % str(resource))

    if 'n_cores_per_node' not in slurm_config[slurm_config['resource']]:
        raise Exception('n_cores_per_node must be set in slurm_config')

    n_cores_per_node = slurm_config[slurm_config['resource']]['n_cores_per_node']
    if type(n_cores_per_node) != int:
        raise Exception('n_cores_per_node (%s) must be an integer.' % str(n_cores_per_node))

    if 'n_nodes_max' not in slurm_config[slurm_config['resource']]:
        raise Exception('n_nodes_max must be set in slurm_config')

    n_nodes_max = slurm_config[slurm_config['resource']]['n_nodes_max']
    if type(n_nodes_max) != int:
        raise Exception('n_nodes_max (%s) must be an integer.' % str(n_nodes_max))

    if 'walltime_max' not in slurm_config[slurm_config['resource']]:
        raise Exception('walltime_max must be set in slurm_config')

    walltime_max = slurm_config[slurm_config['resource']]['walltime_max']
    if type(walltime_max) != str:
        raise Exception('walltime_max (%s) must be formatted as HH:MM:SS' % str(walltime_max))
    if walltime_max.count(':') != 2:
        raise Exception('walltime_max (%s) must be formatted as HH:MM:SS' % str(walltime_max))
    if len(walltime_max) != 8:
        raise Exception('walltime_max (%s) must be formatted as HH:MM:SS' % str(walltime_max))
    for num in walltime_max.split(':'):
        if not num.isdigit():
            raise Exception('walltime_max (%s) must be formatted as HH:MM:SS' % str(walltime_max))

    if type(walltime) != str:
        raise Exception('walltime (%s) must be formatted as HH:MM:SS' % str(walltime_max))
    if walltime.count(':') != 2:
        raise Exception('walltime (%s) must be formatted as HH:MM:SS' % str(walltime))
    if len(walltime) != 8:
        raise Exception('walltime (%s) must be formatted as HH:MM:SS' % str(walltime))
    for num in walltime.split(':'):
        if not num.isdigit():
            raise Exception('walltime (%s) must be formatted as HH:MM:SS' % str(walltime))


def generate_slurm_script(slurm_config_filename, popsycle_config_filename,
                          path_run, output_root,
                          longitude, latitude, area,
                          walltime,
                          n_cores_perform_pop_syn = 1,
                          n_cores_calc_events = 1,
                          n_cores_refine_binary_events = 1,
                          multi_proc_refine_binary_events = True,
                          jobname='default',
                          seed=None, overwrite=False, submitFlag=True,
                          returnJobID=False, dependencyJobID=None,
                          skip_galaxia=False, skip_perform_pop_syn=False,
                          skip_calc_events=False, skip_refine_events=False,
                          skip_refine_binary_events=False,
                          verbose = 0):
    """
    Generates (and possibly submits) the slurm script that
    executes the PopSyCLE pipeline

    Parameters
    ----------
    slurm_config_filename : str
        Name of slurm_config.yaml file containing the slurm parameters
        that will be used the generate the slurm script header.

    popsycle_config_filename : str
        Name of popsycle_config.yaml file containing the PopSyCLE parameters
        that will be passed along to the run_on_slurm.py command in the
        slurm script.

    path_run : str
        Directory containing the parameter file and PopSyCLE output files

    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    walltime : str
        Amount of walltime that the script will request from slurm.
        Format: hh:mm:ss

    Optional Parameters
    -------------------
    jobname : str
        The name of the slurm job and run_popsycle execution file.
        If 'default', the format will be:
            <longitude>_<latitude>_<output_root>

    seed : int
        If non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for SPISEA and PopSyCLE.
        Default None.

    overwrite : bool
        If True, overwrites output files. If False, exists the
        function if output files are already on disk.
        Default is False.

    submitFlag : bool
        If True, script will be submitted to the slurm scheduler
        after being written to disk. If False, it will not be submitted.
        Default is True

    returnJobID : bool
        If True and submitFlag is True,
        function will return the SLURM job id after script submission.
        If False or submitFlag is False,
        function will return None.
        Default is False

    dependencyJobID : int
        If non-None and submitFlag is True, submitted job will only run
        after dependencyJobID is completed with no errors.
        Default is None

    n_cores_perform_pop_syn : int
        Number of cores for executing synthetic.perform_pop_syn
        Default is 1.
        
    n_cores_calc_events : int
        Number of cores for executing synthetic.calc_events
        Default is 1.
        
    n_cores_refine_binary_events : int
        Number of cores for executing synthetic.refine_binary_events
        Default is 1.

    multi_proc_refine_binary_events : bool
        Even if n_proc = 1, a pool is still created. If multi_proc = False,
        instead there is just a for-loop to generate and analyze the lightcurves.
        If multi_proc == False, n_proc must = 1.
        Default is True.

    skip_galaxia : bool
        If True, pipeline will not run Galaxia and assume that the
        resulting ebf file is already present.
        Default is False

    skip_perform_pop_syn : bool
        If True, pipeline will not run perform_pop_syn and assume that
        the resulting h5 file is already present.
        Default is False

    skip_calc_events : bool
        If True, pipeline will not run calc_events and assume that the
        resulting events and blends files are already present.
        Default is False

    skip_refine_events : bool
        If True, pipeline will not run refine_events.
        Default is False

    skip_refine_binary_events : bool
        If True, pipeline will not run refine_binary_events.
        If specified multiplicity is None, will be True.
        Default is False
        
    verbose : int
        Level of debugging information to print to stderr. Set to 0 for minimal
        information. Coarse timing at 2 and fine timing at 4.
        Default is 0.

    Output
    ------
    <output_root>.h5 : hdf5 file
        NOTE: This is what _bin_lb_hdf5 returns.
        An hdf5 file with datasets that correspond to the longitude bin edges,
        latitude bin edges, and the compact objects
        and stars sorted into those bins.

    <path_run>/<output_root>_field_config.yaml : yaml file
        Field config file containing the galactic parameters needed to run pipeline

    <path_run>/run_popsycle_<jobname>.sh : yaml file
        SLURM batch script for submitting job with pipeline run

    Returns
    -------
    None or slurm_jobid : str
        If submitFlag is True and returnJobID is True,
        function returns the SLURM job ID after script submission.
        Otherwise, function returns None

    """
    # Check for files
    if not os.path.exists(slurm_config_filename):
        raise Exception('Slurm configuration file {0} does not exist. '
                        'Write out file using '
                        'run.generate_slurm_config_file '
                        'before proceeding.'.format(slurm_config_filename))
    if not os.path.exists(popsycle_config_filename):
        raise Exception('PopSyCLE configuration file {0} does not exist. '
                        'Write out file using '
                        'run.generate_popsycle_config_file '
                        'before proceeding.'.format(popsycle_config_filename))

    # Enforce popsycle_config_filename is an absolute path
    popsycle_config_filename = os.path.abspath(popsycle_config_filename)
    popsycle_config = load_config_file(popsycle_config_filename)

    # Load multiplicity from popsycle_config
    multiplicity = multiplicity_list[popsycle_config['multiplicity']]
    # Additional multiplicity classes may require a different method of instantiation
    # that would require breaking this out into a separate function
    if multiplicity is not None:
        # These arguments ensure a maximum of triples
        multiplicity = multiplicity(CSF_max=2, companion_max=True)
        hdf5_file_comp = '%s_companions.h5' % output_root
    else:
        skip_refine_binary_events = True
        hdf5_file_comp = None

    # Load the slurm configuration file
    slurm_config = load_config_file(slurm_config_filename)
    
    # Create n_cores dict
    n_cores_popsycle = {'n_cores_perform_pop_syn' : n_cores_perform_pop_syn,
                        'n_cores_calc_events' : n_cores_calc_events,
                        'n_cores_refine_binary_events' : n_cores_refine_binary_events}

    # Check pipeline stages for valid inputs
    _check_slurm_config(slurm_config, walltime)
    if not skip_galaxia:
        _check_run_galaxia(output_root=output_root,
                           longitude=longitude,
                           latitude=latitude,
                           area=area,
                           galaxia_galaxy_model_filename=popsycle_config['galaxia_galaxy_model_filename'],
                           seed=seed)
    if not skip_perform_pop_syn:
        if popsycle_config['bin_edges_number'] == 'None':
            popsycle_config['bin_edges_number'] = None
        _check_perform_pop_syn(ebf_file='test.ebf',
                               output_root=output_root,
                               iso_dir=popsycle_config['isochrones_dir'],
                               IFMR=popsycle_config['IFMR'],
                               bin_edges_number=popsycle_config['bin_edges_number'],
                               BH_kick_speed_mean=popsycle_config['BH_kick_speed_mean'],
                               NS_kick_speed_mean=popsycle_config['NS_kick_speed_mean'],
                               additional_photometric_systems=[popsycle_config['photometric_system']],
                               n_proc=n_cores_perform_pop_syn,
                               binning = popsycle_config['binning'],
                               verbose = verbose,
                               overwrite=overwrite,
                               seed=seed,
                               multiplicity=multiplicity)
    if not skip_calc_events:
        _check_calc_events(hdf5_file='test.h5',
                           output_root2=output_root,
                           radius_cut=popsycle_config['radius_cut'],
                           obs_time=popsycle_config['obs_time'],
                           n_obs=popsycle_config['n_obs'],
                           theta_frac=popsycle_config['theta_frac'],
                           blend_rad=popsycle_config['blend_rad'],
                           n_proc=n_cores_calc_events,
                           overwrite=overwrite,
                           hdf5_file_comp=hdf5_file_comp)
    if not skip_refine_events:
        _check_refine_events(input_root='test',
                             filter_name=popsycle_config['filter_name'],
                             photometric_system=popsycle_config['photometric_system'],
                             red_law=popsycle_config['red_law'],
                             overwrite=overwrite,
                             legacy=False,
                             output_file='default',
                             hdf5_file_comp=hdf5_file_comp,
                             seed=seed)
    if not skip_refine_binary_events:
        refined_events_filename = '{0:s}_refined_events_' \
                              '{1:s}_{2:s}_{3:s}.' \
                              'fits'.format(output_root,
                                            popsycle_config['photometric_system'],
                                            popsycle_config['filter_name'],
                                            popsycle_config['red_law'])
        refined_events_comp_filename = refined_events_filename.replace('.fits', '_companions.fits')
        phot_dir = '%s_bin_phot' % output_root
        _check_refine_binary_events(events=refined_events_filename,
                                    companions=refined_events_comp_filename,
                                    photometric_system=popsycle_config['photometric_system'],
                                    filter_name=popsycle_config['filter_name'],
                                    n_proc=n_cores_refine_binary_events,
                                    overwrite=overwrite,
                                    output_file='default', save_phot=True,
                                    phot_dir=phot_dir, multi_proc=multi_proc_refine_binary_events)


    # Make a run directory for the PopSyCLE output
    path_run = os.path.abspath(path_run)
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # Write a field configuration file to disk in path_run
    config = {'longitude': longitude,
              'latitude': latitude,
              'area': area}
    field_config_filename = '{0}/{1}_field_config.yaml'.format(path_run,
                                                               output_root)
    generate_config_file(field_config_filename, config)

    # Create a slurm jobname base that all stages will be appended to
    if jobname == 'default':
        jobname = 'l%.3f_b%.3f_%s' % (longitude, latitude, output_root)

    ## Bring the slurm_config values into the namespace so that down before
    ## the **locals() command can be executed

    # Path to the python executable
    path_python = slurm_config['path_python']
    # Project account name to charge
    account = slurm_config['account']
    # Queue
    queue = slurm_config['queue']
    # Name of the resource that will be used for the run
    resource = slurm_config['resource']
    # Memory that will be used for the run, in GB
    memory = slurm_config['memory']
    # Maximum number of cores per node
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_nodes_max = slurm_config[resource]['n_nodes_max']
    # Maximum memory on each compute node
    memory_max = slurm_config[resource]['memory_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']
    # Get filepath of the run_on_slurm file
    run_filepath = os.path.dirname(inspect.getfile(load_config_file))
    
    if max(n_cores_popsycle.values()) > n_cores_per_node:
        raise Exception(max(n_cores_popsycle, key = n_cores_popsycle.get) + ' (%s) '
                        'must be less than or equal to '
                        'n_cores_per_node (%s)' % (max(n_cores_popsycle.values()),
                                                   n_cores_per_node))

    if memory > memory_max:
        raise Exception('memory (%s) '
                        'must be less than or equal to '
                        'memory_max (%s)' % (memory, memory_max))

    walltime_s = int(walltime.split(':')[0]) * 3600
    walltime_s += int(walltime.split(':')[1]) * 60
    walltime_s += int(walltime.split(':')[2])
    walltime_max_s = int(walltime_max.split(':')[0]) * 3600
    walltime_max_s += int(walltime_max.split(':')[1]) * 60
    walltime_max_s += int(walltime_max.split(':')[2])
    if walltime_s > walltime_max_s:
        raise Exception('walltime (%s) '
                        'must be less than or equal to '
                        'walltime_max (%s)' % (walltime, walltime_max))

    # Template for writing slurm script. Text must be left adjusted.
    slurm_template = """#!/bin/sh
# Job name
#SBATCH --account={account}
#SBATCH --qos={queue}
#SBATCH --constraint={resource}"""
    if memory < memory_max:
        slurm_template += """
#SBATCH --mem={memory}GB"""
    slurm_template += """
#SBATCH --nodes=1
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
#SBATCH --output={jobname}-%j.out
echo "---------------------------"
echo Longitude = {longitude}
echo Latitude = {latitude}
echo Area = {area}
echo path_run = {path_run}
echo jobname = {jobname} 
echo "Job id = $SLURM_JOBID"
echo "Proc id = $SLURM_PROCID"
hostname
date
echo "---------------------------"

cd {path_run}

"""
    for line in slurm_config['additional_lines']:
        slurm_template += '%s\n' % line
    slurm_template += """
srun -N 1 -n 1 {path_python} {run_filepath}/run.py --output-root={output_root} --field-config-filename={field_config_filename} --popsycle-config-filename={popsycle_config_filename} {optional_cmds}
exitcode=$?
 
date
echo "All done!"
exit $exitcode
"""

    optional_cmds = ''

    # Pass along optional parameters if present
    if not skip_perform_pop_syn:
        optional_cmds += '--n-cores-perform-pop-syn={} '.format(n_cores_perform_pop_syn)
    
    if not skip_calc_events:
        optional_cmds += '--n-cores-calc-events={} '.format(n_cores_calc_events)
    
    if not skip_refine_binary_events:
        optional_cmds += '--n-cores-refine-binary-events={} '.format(n_cores_refine_binary_events)
        optional_cmds += '--multi-proc-refine-binary-events={} '.format(multi_proc_refine_binary_events)

    if overwrite:
        optional_cmds += '--overwrite '

    if seed is not None:
        optional_cmds += '--seed=%i ' % seed

    if skip_galaxia:
        optional_cmds += '--skip-galaxia '

    if skip_perform_pop_syn:
        optional_cmds += '--skip-perform-pop-syn '

    if skip_calc_events:
        optional_cmds += '--skip-calc-events '

    if skip_refine_events:
        optional_cmds += '--skip-refine-events '

    if skip_refine_binary_events:
        optional_cmds += '--skip-refine-binary-events '

    # Populate the mpi_template specified inputs
    job_script = slurm_template.format(**locals())

    # Write the script to the path_run folder
    script_filename = path_run + '/run_popsycle_%s.sh' % (jobname)
    with open(script_filename, 'w') as f:
        f.write(job_script)

    print(f'  {script_filename} written')

    # Submit the job to disk
    slurm_jobid = None
    if submitFlag:
        cwd = os.getcwd()
        os.chdir(path_run)
        if dependencyJobID is not None:
            results = utils.execute('sbatch --dependency=afterok:{0} '
                                    '{1}'.format(dependencyJobID,
                                                 script_filename))
            print('Submitted job {0} to {1} for '
                  '{2} time with {3} dependency'.format(script_filename,
                                                        resource,
                                                        walltime,
                                                        dependencyJobID))
        else:
            results = utils.execute('sbatch {0}'.format(script_filename))
            print('Submitted job {0} to {1} for {2} time'.format(script_filename,
                                                                 resource,
                                                                 walltime))
        stdout, stderr = results
        print('---- Standard Out')
        print(stdout)
        print('---- Standard Err')
        print(stderr)
        print('')
        os.chdir(cwd)

        try:
            slurm_jobid = int(stdout.replace('\n','').split('job')[1])
        except Exception as e:
            slurm_jobid = None

    return slurm_jobid


def tar_run_results(extension_list=['ebf', 'fits', 'h5', 'log', 'out', 'sh', 'txt', 'yaml'],
                    include_bin_phot=True,
                    output_prefix=None):
    """
    Creates a tarball of all results from executing run.py. This script
    assumes that all instances of `generat_slurm_script` were run with
    `path_run` set to subfolders to a common folder, and that this
    script is being run in that common folder. Failure to do so
    will create some funky tarball with files you don't want in it.

    Parameters
    ----------
    extension_list : list
        List of extensions to include in tarball.
        Defaults to all possible extensions.

    include_bin_phot : bool
        If True and running with multiplicity, include lightcurves
        generated by `refine_binary_events`. Default True.

    output_prefix : str
        If not None, add `output_prefix` before name of tarball:
            {output_prefix}_runs.tar

    Output
    ------
    <output_prefix>_runs.tar : tarball file

    Returns
    -------
    None
    """
    extensions = ' '.join(extension_list)
    print(f'Executing "tar_run_results" in CURRENT directory '
          f'for extensions: {extensions}')
    folders = glob.glob('*')
    folders.sort()
    print('-- %i folders gathered' % len(folders))

    fis = []
    for folder in folders:
        for ext in extension_list:
            fis += glob.glob(f'{folder}/*{ext}')
        if include_bin_phot:
            fis += glob.glob(f'{folder}/*bin_phot*/*')


    tar_files_fname = 'tar_files.txt'
    with open(tar_files_fname, 'w') as f:
        for fi in fis:
            f.write('%s\n' % fi)

    print('-- %i files gathered' % len(fis))

    if output_prefix is not None:
        output_fname = f'{output_prefix}_runs.tar'
    else:
        output_fname = 'runs.tar'

    cmd = f'tar -cvf {output_fname} -T {tar_files_fname}'
    print('-- executing tarball creation')
    stdout, stderr = utils.execute(cmd)
    print('-- STDOUT --')
    print(stdout)
    print('-- STDERR --')
    print(stderr)


    os.remove(tar_files_fname)


def run(output_root='root0',
        field_config_filename='field_config.yaml',
        popsycle_config_filename='popsycle_config.yaml',
        n_cores_perform_pop_syn=1,
        n_cores_calc_events=1,
        n_cores_refine_binary_events=1,
        multi_proc_refine_binary_events=True,
        verbose=0,
        seed=None,
        overwrite=False,
        skip_galaxia=False,
        skip_perform_pop_syn=False,
        skip_calc_events=False,
        skip_refine_events=False,
        skip_refine_binary_events=False):

    t0 = time.time()

    # Check for field config file. Exit if not present.
    if not os.path.exists(field_config_filename):
        print("""Error: Field configuration file {0} missing, 
        cannot continue. In order to execute run.py, generate a 
        field configuration file using 
        popsycle.run.generate_field_config_file. 
        Exiting...""".format(field_config_filename))
        sys.exit(1)

    # Check for popsycle config file. Exit if not present.
    if not os.path.exists(popsycle_config_filename):
        print("""Error: popsycle configuration file {0} missing, 
        cannot continue. In order to execute run.py, generate a 
        popsycle configuration file using 
        popsycle.run.generate_popsycle_config_file. 
        Exiting...""".format(popsycle_config_filename))
        sys.exit(1)

    # Load the config files for field parameters
    field_config = load_config_file(field_config_filename)

    # Load the config files for popsycle parameters. If the `bin_edges_number`
    # has been set to the string `None`, instead set it to the boolean None.
    popsycle_config = load_config_file(popsycle_config_filename)
    if popsycle_config['bin_edges_number'] == 'None':
        popsycle_config['bin_edges_number'] = None
    if popsycle_config['photometric_system'] not in synthetic.photometric_system_dict:
        print("""Error: 'photometric_system' in {0} must be a valid option 
        in 'photometric_system_dict'. 
        Exiting...""".format(popsycle_config_filename))
        sys.exit(1)

    # Return the dictionary containing PopSyCLE output filenames
    filename_dict = _return_filename_dict(output_root)

    # Prepare additional_photometric_systems
    additional_photometric_systems = None
    if popsycle_config['photometric_system'] != 'ubv':
        additional_photometric_systems = [popsycle_config['photometric_system']]

    # Load multiplicity from popsycle_config
    multiplicity = multiplicity_list[popsycle_config['multiplicity']]
    # Additional multiplicity classes may require a different method of instantiation
    # that would require breaking this out into a separate function
    if multiplicity is not None:
        # These arguments ensure a maximum of triples
        multiplicity = multiplicity(CSF_max=2, companion_max=True)
        hdf5_file_comp = '%s_companions.h5' % output_root
    else:
        hdf5_file_comp = None

    # Check pipeline stages for valid inputs
    if not skip_galaxia:
        _check_run_galaxia(output_root=output_root,
                           longitude=field_config['longitude'],
                           latitude=field_config['latitude'],
                           area=field_config['area'],
                           galaxia_galaxy_model_filename=popsycle_config['galaxia_galaxy_model_filename'],
                           seed=seed)
    if not skip_perform_pop_syn:
        _check_perform_pop_syn(ebf_file=filename_dict['ebf_filename'],
                               output_root=output_root,
                               iso_dir=popsycle_config['isochrones_dir'],
                               IFMR=popsycle_config['IFMR'],
                               bin_edges_number=popsycle_config['bin_edges_number'],
                               BH_kick_speed_mean=popsycle_config['BH_kick_speed_mean'],
                               NS_kick_speed_mean=popsycle_config['NS_kick_speed_mean'],
                               additional_photometric_systems=additional_photometric_systems,
                               verbose=verbose,
                               n_proc=n_cores_perform_pop_syn,
                               binning=popsycle_config['binning'],
                               overwrite=overwrite,
                               seed=seed,
                               multiplicity=multiplicity)
    if not skip_calc_events:
        _check_calc_events(hdf5_file=filename_dict['hdf5_filename'],
                           output_root2=output_root,
                           radius_cut=popsycle_config['radius_cut'],
                           obs_time=popsycle_config['obs_time'],
                           n_obs=popsycle_config['n_obs'],
                           theta_frac=popsycle_config['theta_frac'],
                           blend_rad=popsycle_config['blend_rad'],
                           n_proc=n_cores_calc_events,
                           overwrite=overwrite,
                           hdf5_file_comp=hdf5_file_comp)
    if not skip_refine_events:
        _check_refine_events(input_root=output_root,
                             filter_name=popsycle_config['filter_name'],
                             photometric_system=popsycle_config['photometric_system'],
                             red_law=popsycle_config['red_law'],
                             overwrite=overwrite,
                             legacy=False,
                             output_file='default',
                             hdf5_file_comp=hdf5_file_comp,
                             seed=seed)
    if not skip_refine_binary_events:
        refined_events_filename = '{0:s}_refined_events_' \
                              '{1:s}_{2:s}_{3:s}.' \
                              'fits'.format(output_root,
                                            popsycle_config['photometric_system'],
                                            popsycle_config['filter_name'],
                                            popsycle_config['red_law'])
        refined_events_comp_filename = refined_events_filename.replace('.fits', '_companions.fits')
        phot_dir = '%s_bin_phot' % output_root
        _check_refine_binary_events(events=refined_events_filename,
                                    companions=refined_events_comp_filename,
                                    photometric_system=popsycle_config['photometric_system'],
                                    filter_name=popsycle_config['filter_name'],
                                    n_proc=n_cores_refine_binary_events,
                                    overwrite=overwrite,
                                    output_file='default', save_phot=True,
                                    phot_dir=phot_dir, multi_proc=multi_proc_refine_binary_events)

    if not skip_galaxia:
        # Remove Galaxia output if already exists and overwrite=True
        if _check_for_output(filename_dict['ebf_filename'],
                             overwrite):
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)

        # Run Galaxia
        print('-- Running Galaxia')
        synthetic.run_galaxia(output_root=output_root,
                              longitude=field_config['longitude'],
                              latitude=field_config['latitude'],
                              area=field_config['area'],
                              galaxia_galaxy_model_filename=popsycle_config['galaxia_galaxy_model_filename'],
                              seed=seed)

    if not skip_perform_pop_syn:
        # Remove perform_pop_syn output if already exists and overwrite=True
        if _check_for_output(filename_dict['hdf5_filename'],
                             overwrite):
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)

        # Run perform_pop_syn
        print('-- Executing perform_pop_syn')
        synthetic.perform_pop_syn(
            ebf_file=filename_dict['ebf_filename'],
            output_root=output_root,
            iso_dir=popsycle_config['isochrones_dir'],
            IFMR=popsycle_config['IFMR'],
            bin_edges_number=popsycle_config['bin_edges_number'],
            BH_kick_speed_mean=popsycle_config['BH_kick_speed_mean'],
            NS_kick_speed_mean=popsycle_config['NS_kick_speed_mean'],
            additional_photometric_systems=additional_photometric_systems,
            n_proc=n_cores_perform_pop_syn,
            overwrite=overwrite,
            seed=seed,
            multiplicity=multiplicity)

    if not skip_calc_events:
        # Remove calc_events output if already exists and overwrite=True
        if _check_for_output(filename_dict['events_filename'],
                             overwrite):
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)
        if _check_for_output(filename_dict['blends_filename'],
                             overwrite):
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)

        # Run calc_events
        print('-- Executing calc_events')
        synthetic.calc_events(hdf5_file=filename_dict['hdf5_filename'],
                              output_root2=output_root,
                              radius_cut=popsycle_config['radius_cut'],
                              obs_time=popsycle_config['obs_time'],
                              n_obs=popsycle_config['n_obs'],
                              theta_frac=popsycle_config['theta_frac'],
                              blend_rad=popsycle_config['blend_rad'],
                              n_proc=n_cores_calc_events,
                              overwrite=overwrite,
                              hdf5_file_comp=hdf5_file_comp)

        # Write a fle to disk stating that there are no events if
        # calc_events does not produce an events file
        if not os.path.exists(filename_dict['events_filename']):
            Path(filename_dict['noevents_filename']).touch()
            print('No events present, skipping refine_events')
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(0)

    refined_events_filename = '{0:s}_refined_events_' \
                              '{1:s}_{2:s}_{3:s}.' \
                              'fits'.format(output_root,
                                            popsycle_config['photometric_system'],
                                            popsycle_config['filter_name'],
                                            popsycle_config['red_law'])
    if not skip_refine_events:
        # Remove refine_events output if already exists and overwrite=True
        if _check_for_output(refined_events_filename, overwrite):
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)

        # Run refine_events
        print('-- Executing refine_events')
        synthetic.refine_events(input_root=output_root,
                                filter_name=popsycle_config['filter_name'],
                                photometric_system=popsycle_config['photometric_system'],
                                red_law=popsycle_config['red_law'],
                                overwrite=overwrite,
                                output_file='default',
                                hdf5_file_comp=hdf5_file_comp,
                                seed=seed)

    if multiplicity is not None and not skip_refine_binary_events:
        if not os.path.exists(refined_events_filename):
            print('Refined events %s missing and therefore '
                  'cannot run refine_binary_events. Skipping...'
                  % refined_events_filename)
            t1 = time.time()
            print('run.py runtime : {0:f} s'.format(t1 - t0))
            sys.exit(1)

        refined_events_comp_filename = refined_events_filename.replace('.fits', '_companions.fits')
        phot_dir = '%s_bin_phot' % output_root
        synthetic.refine_binary_events(events=refined_events_filename,
                                       companions=refined_events_comp_filename,
                                       photometric_system=popsycle_config['photometric_system'],
                                       filter_name=popsycle_config['filter_name'],
                                       n_proc=n_cores_refine_binary_events,
                                       overwrite=overwrite,
                                       output_file='default', save_phot=True,
                                       phot_dir=phot_dir,
                                       multi_proc=multi_proc_refine_binary_events)

    t1 = time.time()
    print('run.py runtime : {0:f} s'.format(t1 - t0))


def main():
    description_str = """
    Run the PopSyCLE pipeline. This executable can be either 
    run by slurm scripts generated by `generate_slurm_scripts` or from the 
    command line.

    Script must be executed in a folder containing a field_config file and 
    point to a popsycle_config file both generated by 
    `popsycle.slurm.generate_config_file`.
    """

    parser = argparse.ArgumentParser(description=description_str,
                                     formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group(title='Required')
    required.add_argument('--output-root', type=str,
                          help='Base filename of the output files. '
                               'Default: root0',
                          default='root0')
    required.add_argument('--field-config-filename', type=str,
                          help='Name of configuration file containing '
                               'the field parameters. '
                               'Default: field_config.yaml',
                          default='field_config.yaml')
    required.add_argument('--popsycle-config-filename', type=str,
                          help='Name of configuration file containing '
                               'the PopSyCLE parameters. '
                               'Default: popsycle_config.yaml',
                          default='popsycle_config.yaml')

    optional = parser.add_argument_group(title='Optional')
    optional.add_argument('--n-cores-perform-pop-syn', type=int,
                          help='Number of cores to use in the perform_pop_syn '
                               'function.'
                               'Default is --n-cores=1 or serial processing.',
                          default=1)
    optional.add_argument('--n-cores-calc-events', type=int,
                          help='Number of cores to use in the calc_events '
                               'function. '
                               'Default is --n-cores=1 or serial processing.',
                          default=1)
    optional.add_argument('--n-cores-refine-binary-events', type=int,
                          help='Number of cores to use in the refine_binary_events '
                               'function. '
                               'Default is --n-cores=1 or serial processing.',
                          default=1)
    optional.add_argument('--multi-proc-refine-binary-events', type=bool,
                          help='Controls multi processing for refine bianry events '
                          'even if n-cores=1',
                          default=True)
    optional.add_argument('--seed', type=int,
                          help='Set a seed for all PopSyCLE functions with '
                               'randomness, which are running Galaxia and '
                               'SPISEA. Setting this flag guarantees '
                               'identical output and is useful for debugging.',
                          default=None)
    optional.add_argument('--overwrite',
                          help="Overwrite all output files.",
                          action='store_true')
    optional.add_argument('--skip-galaxia',
                          help="Skip running galaxia.",
                          action='store_true')
    optional.add_argument('--skip-perform-pop-syn',
                          help="Skip running perform_pop_syn.",
                          action='store_true')
    optional.add_argument('--skip-calc-events',
                          help="Skip running calc_events.",
                          action='store_true')
    optional.add_argument('--skip-refine-events',
                          help="Skip running refine_events.",
                          action='store_true')
    optional.add_argument('--skip-refine-binary-events',
                          help="Skip running refine_binary_events.",
                          action='store_true')
    args = parser.parse_args()

    run(output_root=args.output_root,
        field_config_filename=args.field_config_filename,
        popsycle_config_filename=args.popsycle_config_filename,
        n_cores_perform_pop_syn=args.n_cores_perform_pop_syn,
        n_cores_calc_events=args.n_cores_calc_events,
        n_cores_refine_binary_events=args.n_cores_refine_binary_events,
        multi_proc_refine_binary_events=args.multi_proc_refine_binary_events,
        seed=args.seed,
        overwrite=args.overwrite,
        skip_galaxia=args.skip_galaxia,
        skip_perform_pop_syn=args.skip_perform_pop_syn,
        skip_calc_events=args.skip_calc_events,
        skip_refine_events=args.skip_refine_events,
        skip_refine_binary_events=args.skip_refine_binary_events)


if __name__ == '__main__':
    main()
