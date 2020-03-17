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
from popsycle import synthetic
from popsycle import utils


def _return_filename_dict(output_root, add_pbh_flag=False):
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

    # Append _pbh to filenames if add_pbh will be run
    if add_pbh_flag:
        events_filename = '%s_pbh_events.fits' % output_root
        blends_filename = '%s_pbh_blends.fits' % output_root
        noevents_filename = '%s_pbh_NOEVENTS.txt' % output_root
    else:
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


def load_config(config_filename):
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
        config = yaml.safe_load(f)
    return config


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


def generate_slurm_config_file(path_python, account, queue,
                               resource, n_cores_per_node, n_nodes_max,
                               walltime_max, additional_lines,
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

    n_cores_per_node : int
        Number of cores in each node of the compute resource

    n_nodes_max : int
        Total number of nodes in the compute resource

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
              'additional_lines': additional_lines,
              resource: {'n_cores_per_node': n_cores_per_node,
                         'n_nodes_max': n_nodes_max,
                         'walltime_max': walltime_max}}
    generate_config_file(config_filename, config)


def generate_popsycle_config_file(radius_cut, obs_time,
                                  n_obs, theta_frac, blend_rad,
                                  isochrones_dir,
                                  bin_edges_number,
                                  BH_kick_speed_mean, NS_kick_speed_mean,
                                  photometric_system,
                                  filter_name, red_law,
                                  config_filename='popsycle_config.yaml'):
    """
    Save popsycle configuration parameters from a dictionary into a yaml file

    Parameters
    ----------
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

    isochrones_dir : str
        Directory for PyPopStar isochrones

    bin_edges_number : int
        Number of edges for the bins (bins = bin_edges_number - 1)
        Total number of bins is (bin_edges_number - 1)**2

    BH_kick_speed_mean : float
        Mean of the birth kick speed of BH (in km/s) maxwellian distrubution.

    NS_kick_speed_mean : float
        Mean of the birth kick speed of NS (in km/s) maxwellian distrubution.

    photometric_system : str
        The name of the photometric system in which the filter exists.

    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    red_law : str
        The name of the reddening law to use from PopStar.

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

    config = {'radius_cut': radius_cut,
              'obs_time': obs_time,
              'n_obs': n_obs,
              'theta_frac': theta_frac,
              'blend_rad': blend_rad,
              'isochrones_dir': isochrones_dir,
              'bin_edges_number': bin_edges_number,
              'BH_kick_speed_mean': BH_kick_speed_mean,
              'NS_kick_speed_mean': NS_kick_speed_mean,
              'photometric_system': photometric_system,
              'filter_name': filter_name,
              'red_law': red_law}
    generate_config_file(config_filename, config)


def generate_pbh_config_file(config_filename, fdm, pbh_mass, r_max, r_s, gamma, v_esc, rho_0, n_lin):
    """
    Save PBH configuration parameters into a yaml file

    Parameters
    ----------

    fdm : float
        Fraction of dark matter.
        The fraction of dark matter that you want to consist of PBHs.

    pbh_mass : int
        The single mass that all PBHs will have (in units of Msun).

    r_max : float
        The maximum radius (in kpc) from the Earth that you want to find PBHs.

    r_s: float
        The scale radius of the Milky Way (in units of kpc). r_s = r_vir / c (virial radius / concentration index) 

    gamma: float
        The inner slope of the MW dark matter halo as described in LaCroix 2018.
        Gamma goes into the determination of the velocities and each value returns a slightly different distribution.

    v_esc: int
        The escape velocity of the Milky Way in km/s.
        v_esc is used in calculating the velocities.

    rho_0: float
        The initial density that will be used in the NFW profile equations (in units of Msun/pc^3).

    n_lin: int
        The number of times you want the density determined along the line of sight when calculating PBH positions
        Defaults to 1000. Will need to make large if you are closer to the galactic center.

    Output
    ------
    None
    """

    config = {'fdm': fdm,
              'pbh_mass': pbh_mass,
              'r_max': r_max,
              'r_s': r_s,
              'gamma': gamma,
              'v_esc': v_esc,
              'rho_0':rho_0,
              'n_lin':n_lin}
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


def generate_slurm_script(slurm_config_filename, popsycle_config_filename,
                          path_run, output_root,
                          longitude, latitude, area,
                          n_cores_calc_events,
                          walltime,
                          seed=None, overwrite=False, submitFlag=True,
                          skip_galaxia=False, skip_perform_pop_syn=False,
                          skip_calc_events=False, skip_refine_events=False,
                          pbh_config_filename=None):
    """
    Generates the slurm script that executes the PopSyCLE pipeline

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

    n_cores_calc_events : int
        Number of cores for executing synthetic.calc_events

    walltime : str
        Amount of walltime that the script will request from slurm.
        Format: hh:mm:ss

    Optional Parameters
    -------------------
    pbh_config_filename : str
        Name of pbh_config.yaml file containing the PBH parameters
        that will be passed along to the run_on_slurm.py command in the
        slurm script.

    seed : int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output for PyPopStar and PopSyCLE.
        Default None.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    submitFlag : bool
        If set to True, script will be submitted to the slurm scheduler
        after being written to disk. If set to False, it will not be submitted.
        Default is True

    skip_galaxia : bool
        If set to True, pipeline will not run Galaxia and assume that the
        resulting ebf file is already present.
        Default is False

    skip_perform_pop_syn : bool
        If set to True, pipeline will not run perform_pop_syn and assume that
        the resulting h5 file is already present.
        Default is False

    skip_calc_events : bool
        If set to True, pipeline will not run calc_events and assume that the
        resulting events and blends files are already present.
        Default is False

    skip_refine_events : bool
        If set to True, pipeline will not run refine_events.
        Default is False

    Output
    ------
    None

    """
    # Check for files
    if not os.path.exists(slurm_config_filename):
        raise Exception('Slurm configuration file {0} does not exist. '
                        'Write out file using synthetic.generate_config_file '
                        'before proceeding.'.format(slurm_config_filename))
    if not os.path.exists(popsycle_config_filename):
        raise Exception('PopSyCLE configuration file {0} does not exist. '
                        'Write out file using synthetic.generate_config_file '
                        'before proceeding.'.format(popsycle_config_filename))

    # Enforce popsycle_config_filename is an absolute path
    popsycle_config_filename = os.path.abspath(popsycle_config_filename)

    # Make a run directory for the PopSyCLE output
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # Write a field configuration file to disk in path_run
    config = {'longitude': longitude,
              'latitude': latitude,
              'area': area}
    field_config_filename = '{0}/field_config.{1}.yaml'.format(path_run,
                                                               output_root)
    generate_config_file(field_config_filename, config)

    # Load the slurm configuration file
    slurm_config = load_config(slurm_config_filename)

    # Create a slurm jobname base that all stages will be appended to
    jobname = 'l%.1f_b%.1f_%s' % (longitude, latitude, output_root)

    ## Bring the slurm_config values into the namespace so that down before
    ## the **locals() command can be executed

    # Path to the python executable
    path_python = slurm_config['path_python']
    # Project account name to charge
    account = slurm_config['account']
    # Queue
    queue = slurm_config['queue']
    # Name of the resource that will be ussed for the run
    resource = slurm_config['resource']
    # Maximum number of ores per node
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_nodes_max = slurm_config[resource]['n_nodes_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']
    # Get filepath of the run_on_slurm file
    run_filepath = os.path.dirname(inspect.getfile(load_config))

    # Template for writing slurm script. Text must be left adjusted.
    slurm_template = """#!/bin/sh
# Job name
#SBATCH --account={account}
#SBATCH --qos={queue}
#SBATCH --constraint={resource}
#SBATCH --nodes=1
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
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
"""
    for line in slurm_config['additional_lines']:
        slurm_template += '%s\n' % line
    slurm_template += """
cd {path_run}
srun -N 1 -n 1 {path_python} {run_filepath}/run.py --output-root={output_root} --field-config-filename={field_config_filename} --popsycle-config-filename={popsycle_config_filename} --n-cores-calc-events={n_cores_calc_events} {optional_cmds} 
date
echo "All done!"
"""

    # Check that the specified number of cores does not exceed the resource max
    if n_cores_calc_events > n_cores_per_node:
        print('Error: specified number of cores exceeds limit. Exiting...')
        return None

    optional_cmds = ''

    # Pass along optional parameters if present
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

    if pbh_config_filename:
        optional_cmds += '--pbh-config-filename={0}'.format(pbh_config_filename)

    # Populate the mpi_template specified inputs
    job_script = slurm_template.format(**locals())

    # Write the script to the path_run folder
    script_filename = path_run + '/run_popsycle_%s.sh' % (jobname)
    with open(script_filename, 'w') as f:
        f.write(job_script)

    # Submit the job to disk
    if submitFlag:
        os.chdir(path_run)
        stdout, stderr = utils.execute('sbatch {0}'.format(script_filename))
        print('Submitted job {0} to {1} for {2} time'.format(script_filename,
                                                             resource,
                                                             walltime))
        print('---- Standard Out')
        print(stdout)
        print('---- Standard Err')
        print(stderr)
        print('')


def run():
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
    required.add_argument('--n-cores-calc-events', type=int,
                          help='Number of cores to use in the calc_events '
                               'function (the only piece of the '
                               'PopSyCLE pipeline that uses multiprocessing). '
                               'Default is --n-cores=1 or serial processing.',
                          default=1)

    optional = parser.add_argument_group(title='Optional')
    optional.add_argument('--seed', type=int,
                          help='Set a seed for all PopSyCLE functions with '
                               'randomness, including running Galaxia and '
                               'PyPopStar. Setting this flag guarantees '
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
    optional.add_argument('--pbh-config-filename', type=str,
                          help='Name of configuration file containing '
                               'pbh inputs. Default if needed is: '
                               'pbh_config.yaml')

    args = parser.parse_args()


    # Check for field config file. Exit if not present.
    if not os.path.exists(args.field_config_filename):
        print("""Error: Field configuration file {0} missing, 
        cannot continue. In order to execute run.py, generate a 
        field configuration file using 
        popsycle.synthetic.generate_field_config_file. 
        Exiting...""".format(args.field_config_filename))
        sys.exit(1)

    # Check for popsycle config file. Exit if not present.
    if not os.path.exists(args.popsycle_config_filename):
        print("""Error: popsycle configuration file {0} missing, 
        cannot continue. In order to execute run.py, generate a 
        popsycle configuration file using 
        popsycle.synthetic.generate_popsycle_config_file. 
        Exiting...""".format(args.popsycle_config_filename))
        sys.exit(1)

    # Load the config files for field parameters
    field_config = load_config(args.field_config_filename)

    # Load the config files for popsycle parameters. If the `bin_edges_number`
    # has been set to the string `None`, instead set it to the boolean None.
    popsycle_config = load_config(args.popsycle_config_filename)
    if popsycle_config['bin_edges_number'] == 'None':
        popsycle_config['bin_edges_number'] = None
    if popsycle_config['photometric_system'] not in synthetic.photometric_system_dict:
        print("""Error: 'photometric_system' in {0} must be a valid option 
        in 'photometric_system_dict'. 
        Exiting...""".format(args.popsycle_config_filename))
        sys.exit(1)

    # Return the dictionary containing PopSyCLE output filenames
    # Append '_pbh' to filenames if add_pbh will be run
    if args.pbh_config_filename is not None:
        add_pbh_flag = True
    else:
        add_pbh_flag = False
    filename_dict = _return_filename_dict(args.output_root,
                                          add_pbh_flag=add_pbh_flag)

    # Prepare additional_photometric_systems
    additional_photometric_systems = None
    if popsycle_config['photometric_system'] != 'ubv':
        additional_photometric_systems = [popsycle_config['photometric_system']]

    if not args.skip_galaxia:
        # Remove Galaxia output if already exists and overwrite=True
        if _check_for_output(filename_dict['ebf_filename'],
                             args.overwrite):
            sys.exit(1)

        # Run Galaxia
        print('-- Running Galaxia')
        synthetic.run_galaxia(output_root=args.output_root,
                              longitude=field_config['longitude'],
                              latitude=field_config['latitude'],
                              area=field_config['area'],
                              seed=args.seed)

    if not args.skip_perform_pop_syn:
        # Remove perform_pop_syn output if already exists and overwrite=True
        if _check_for_output(filename_dict['hdf5_filename'],
                             args.overwrite):
            sys.exit(1)

        # Run perform_pop_syn
        print('-- Executing perform_pop_syn')
        synthetic.perform_pop_syn(
            ebf_file=filename_dict['ebf_filename'],
            output_root=args.output_root,
            iso_dir=popsycle_config['isochrones_dir'],
            bin_edges_number=popsycle_config['bin_edges_number'],
            BH_kick_speed_mean=popsycle_config['BH_kick_speed_mean'],
            NS_kick_speed_mean=popsycle_config['NS_kick_speed_mean'],
            additional_photometric_systems=additional_photometric_systems,
            overwrite=args.overwrite,
            seed=args.seed)

    # Append '_pbh' to output_root if add_pbh will be run
    if args.pbh_config_filename is not None:
        output_root = '%s_pbh' % args.output_root
    else:
        output_root = args.output_root

    # only do stuff if optional config file for pbhs was provided
    if args.pbh_config_filename is not None:
        # Check for pbh config file. Exit if not present.
        if not os.path.exists(args.pbh_config_filename):
            print("""Error: PBH configuration file {0} missing, 
            cannot continue. In order to execute run.py, generate a 
            PBH configuration file using 
            popsycle.synthetic.generate_pbh_config_file. 
            Exiting...""".format(args.pbh_config_filename))
            sys.exit(1)

        pbh_config = synthetic.load_config(args.pbh_config_filename)

        # Check if .h5 file exists from perform popsyn, use as input for following function
        if not os.path.exists(filename_dict['hdf5_filename']):
          print("""Error: H5 file was not created properly by 
            synthetic.perform_pop_syn""")
          sys.exit(1)

        synthetic.add_pbh(hdf5_file=filename_dict['hdf5_filename'],
                          ebf_file=filename_dict['ebf_filename'],
                          output_root2=output_root,
                          fdm=pbh_config['fdm'],
                          pbh_mass=pbh_config['pbh_mass'],
                          r_max=pbh_config['r_max'],
                          r_s=pbh_config['r_s'],
                          gamma=pbh_config['gamma'],
                          v_esc=pbh_config['v_esc'],
                          rho_0=pbh_config['rho_0'],
                          n_lin=pbh_config['n_lin'],
                          overwrite=args.overwrite,
                          seed=args.seed)

    if not args.skip_calc_events:
        # Remove calc_events output if already exists and overwrite=True
        if _check_for_output(filename_dict['events_filename'],
                             args.overwrite):
            sys.exit(1)
        if _check_for_output(filename_dict['blends_filename'],
                             args.overwrite):
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
                              n_proc=args.n_cores_calc_events,
                              seed=args.seed,
                              overwrite=args.overwrite)

        # Write a fle to disk stating that there are no events if
        # calc_events does not produce an events file
        if not os.path.exists(filename_dict['events_filename']):
            Path(filename_dict['noevents_filename']).touch()
            print('No events present, skipping refine_events')
            sys.exit(0)

    if not args.skip_refine_events:
        # Remove refine_events output if already exists and overwrite=True
        filename = '{0:s}_refined_events_{1:s}_{2:s}.' \
                   'fits'.format(output_root,
                                 popsycle_config['filter_name'],
                                 popsycle_config['red_law'])
        if _check_for_output(filename, args.overwrite):
            sys.exit(1)

        # Run refine_events
        print('-- Executing refine_events')
        synthetic.refine_events(input_root=output_root,
                                filter_name=popsycle_config['filter_name'],
                                photometric_system=popsycle_config['photometric_system'],
                                red_law=popsycle_config['red_law'],
                                overwrite=args.overwrite,
                                output_file='default')


if __name__ == '__main__':
    run()
