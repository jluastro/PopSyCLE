#! /usr/bin/env python
"""
run_on_slurm.py
Generate and execute scripts for parallel execution of PopSyCLE runs. Scripts
created will be formatted for submission to a SLURM scheduler. This file
also serves as the executed of the PopSyCLE pipeline by those slurm scripts.
"""

import subprocess
import yaml
import os
from pathlib import Path
import argparse
from popsycle import synthetic


def execute(cmd,
            shell=False):
    """Executes a command line instruction, captures the stdout and stderr

    Args:
        cmd : str
            Command line instruction, including any executables and parameters.
        shell : bool
            Determines if the command is run through the shell.

    Returns:
        stdout : str
            Contains the standard output of the executed process.
        stderr : str
            Contains the standard error of the executed process.

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

    return stdout, stderr


def generate_galactic_config_file(path_run, output_root,
                                  longitude, latitude, area):
    """
    Generate a yaml file that contains the microlensing parameters specific to
    a PopSyCLE run.

    Parameters
    ----------
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


    Output
    ------
    None

    """
    # Name the config file using the output root
    config_filename = '{0}/galactic_config.{1}.yaml'.format(path_run,
                                                            output_root)
    # Write galactic parameters to the config file
    with open(config_filename, 'w') as f:
        f.write('longitude: %f\n' % longitude)
        f.write('latitude: %f\n' % latitude)
        f.write('area: %f\n' % area)


def generate_stage_script(stage, slurm_config, popsycle_config_filename,
                          path_run, output_root,
                          jobname_base, walltime,
                          N_nodes=1, N_cores=1,
                          previous_slurm_job_id=None,
                          submitFlag=True, debugFlag=False):
    """
    Generate a slurm script that executes a stage of the PopSyCLE pipeline

    Parameters
    ----------
    stage : int
        Number 1, 2 or 3 indicating the stage of the PopSyCLE pipeline.
        Stage 1: (serial)
            - Galaxia
            - synthetic.perform_pop_syn
        Stage 2: (parallel)
            - synthetic.calc_events
        Stage 3: (serial)
            - synthetic.refine_events

    slurm_config : dict
        Loaded from a slurm_config.yaml file

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

    jobname_base : str
        Name of slurm jobname, appended by the stage of the script.
        Default format is 'l10.00_b1.00' for a position of (l,b) = (10, 1)

    walltime : str
        Amount of walltime that the script will request from slurm.
        Format must be 'hh:mm:ss'

    N_nodes : int
        Number of nodes to run the stage

    N_cores : int
        Number of cores to run the stage

    previous_slurm_job_id : str
        Slurm Job ID of the previous stage which must be completed before
        this stage can be exected

    submitFlag : bool
        If set to True, scripts will be submitted to the slurm scheduler after
        being written to disk. If set to False, they will not be submitted.
        Default is True

    debugFlag : bool
        If set to True, scripts will be run with a fixed seed that produces
        identical output. If set to False, a random seed will be selected.
        Default is False

    Output
    ------
    None

    """
    # Check that the stage is one of the acceptable inputs
    if stage not in [1, 2, 3]:
        print('Error: stage must be one of [1, 2, 3]. Exiting...')
        os.exit(1)

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
    # Maximum number of cores per node
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_nodes_max = slurm_config[resource]['n_nodes_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']
    # Get filepath of the run_on_slurm file
    run_on_slurm_filepath = os.path.abspath(__file__)
    # Create jobname
    jobname = '%s_s%i' % (jobname_base, stage)

    # Template for writing slurm script. Text must be left adjusted.
    mpi_template = """#!/bin/sh
# Job name
#SBATCH --account={account}
#SBATCH --qos={queue}
#SBATCH --constraint={resource}
#SBATCH --nodes={n_nodes}
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
echo "---------------------------"
date
echo "Job id = $SLURM_JOBID"
echo "Proc id = $SLURM_PROCID"
hostname
echo "---------------------------"
module load cray-hdf5/1.10.5.2
export HDF5_USE_FILE_LOCKING=FALSE
cd {path_run}
srun -N{N_nodes} -n{N_cores} {path_python} {run_on_slurm_filepath} --output-root={output_root} --stage={stage} --microlensing-config-filename={popsycle_config_filename} {debug_cmd} 
date
echo
"All done!"
"""

    # Check that the specified number of nodes does not exceed the resource max
    if N_nodes > n_nodes_max:
        print('Error: specified number of nodes exceeds limit. Exiting...')
        os.exit(1)
    # Check that the specified number of cores does not exceed the resource max
    if N_cores > n_cores_per_node:
        print('Error: specified number of cores exceeds limit. Exiting...')
        os.exit(1)

    # Make a run directory for the PopSyCLE output
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # If debugFlag == True, add debug flag to srun command
    if debugFlag:
        debug_cmd = '--debug'
    else:
        debug_cmd = ''

    # Populate the mpi_template specified inputs
    job_script = mpi_template.format(**locals())

    # Write the script to the path_run folder
    if not path_run.endswith('/'):
        path_run = path_run + '/'
    script_filename = path_run + 'run_popsycle_{0}.sh'.format(stage)
    with open(script_filename, 'w') as f:
        f.write(job_script)

    # Submit the job to disk
    if submitFlag:
        # Add a dependency if previous_slurm_id is not none
        if previous_slurm_job_id is None:
            stdout, stderr = execute('sbatch {0}'.format(script_filename))
        else:
            stdout, stderr = execute(
                'sbatch --dependency=afterok:{0} {1}'.format(
                    previous_slurm_job_id, script_filename))

        if debugFlag:
            print('** Standard Out **')
            print(stdout)
            print('** Standard Err **')
            print(stderr)

        print('Submitted job {0} to {1}'.format(script_filename, resource))

        # Extract slurm job ID of submitted job from stdout
        previous_slurm_job_id = stdout.decode().replace('\n', '').split()[-1]
    else:
        previous_slurm_job_id = None

    return previous_slurm_job_id


def generate_slurm_scripts(slurm_config_filename, popsycle_config_filename,
                           path_run, output_root,
                           longitude, latitude, area,
                           N_nodes_calc_events, N_cores_calc_events,
                           walltime_stage1, walltime_stage2,
                           walltime_stage3, submitFlag=True, debugFlag=False):
    """
    Generates all stages of slurm scripts that executes the PopSyCLE pipeline

    Parameters
    ----------
    stage : int
        Number 1, 2 or 3 indicating the stage of the PopSyCLE pipeline.
        Stage 1: (serial)
            - Galaxia
            - synthetic.perform_pop_syn
        Stage 2: (parallel)
            - synthetic.calc_events
        Stage 3: (serial)
            - synthetic.refine_events

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

    N_nodes_calc_events : int
        Number of nodes for stage 2 where synthetic.calc_events is executed

    N_cores_calc_events : int
        Number of cores for stage 2 where synthetic.calc_events is executed

    walltime_stage1 : str
        Amount of walltime that the script will request from slurm for stage 1
        Format must be 'hh:mm:ss'

    walltime_stage2 : str
        Amount of walltime that the script will request from slurm for stage 2
        Format must be 'hh:mm:ss'

    walltime_stage3 : str
        Amount of walltime that the script will request from slurm for stage 3
        Format must be 'hh:mm:ss'

    submitFlag : bool
        If set to True, scripts will be submitted to the slurm scheduler after
        being written to disk. If set to False, they will not be submitted.
        Default is True

    debugFlag : bool
        If set to True, scripts will be run with a fixed seed that produces
        identical output. If set to False, a random seed will be selected.
        Default is False

    Output
    ------
    None

    """
    # Write a galactic configuration file to disk in path_run
    generate_galactic_config_file(path_run, output_root,
                                  longitude, latitude, area)

    # Load the slurm configuration file
    with open(slurm_config_filename, 'r') as f:
        slurm_config = yaml.safe_load(f)

    # Create a slurm jobname base that all stages will be appended to
    jobname_base = 'l%.1f_b%.1f' % (longitude, latitude)

    # Create the first stage script and submit it if submitFlag == True
    slurm_job_id1 = generate_stage_script(1, slurm_config,
                                          popsycle_config_filename,
                                          path_run, output_root,
                                          jobname_base, walltime_stage1,
                                          submitFlag=submitFlag,
                                          debugFlag=debugFlag)
    # Create the second stage script and submit it if submitFlag == True
    slurm_job_id2 = generate_stage_script(2, slurm_config,
                                          popsycle_config_filename,
                                          path_run, output_root,
                                          jobname_base, walltime_stage2,
                                          N_nodes=N_nodes_calc_events,
                                          N_cores=N_cores_calc_events,
                                          previous_slurm_job_id=slurm_job_id1,
                                          submitFlag=submitFlag,
                                          debugFlag=debugFlag)
    # Create the third stage script and submit it if submitFlag == True
    _ = generate_stage_script(3, slurm_config, popsycle_config_filename,
                              path_run, output_root,
                              jobname_base, walltime_stage3,
                              previous_slurm_job_id=slurm_job_id2,
                              submitFlag=submitFlag,
                              debugFlag=debugFlag)


def load_galactic_config(output_root):
    """
    Load the galactic parameters from the yaml file

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
    galactic_config : dict
        Dictionary containing the galactic parameters for the PopSyCLE run

    """
    # Load the configuration file containing galactic parameters
    config_filename = 'galactic_config.{0}.yaml'.format(output_root)
    with open(config_filename, 'r') as f:
        galactic_config = yaml.safe_load(f)
    return galactic_config


def load_popsycle_config(popsycle_config_filename):
    """
    Load the PopSyCLE parameters from the yaml file

    Parameters
    ----------
    popsycle_config_filename : str
        Name of popsycle_config.yaml file containing the PopSyCLE parameters
        that will be passed along to the run_on_slurm.py command in the
        slurm script.

    Output
    ------
    popsycle_config : dict
        Dictionary containing the PopSyCLE parameters for the PopSyCLE run

    """
    # Load the configuration file containing popsycle parameters
    with open(popsycle_config_filename, 'r') as f:
        popsycle_config = yaml.safe_load(f)
    return popsycle_config


def return_filename_dict(output_root):
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


def run_stage1(output_root,
               galactic_config, popsycle_config,
               filename_dict, debugFlag=False):
    """
    Run stage 1 of the PopSyCLE pipeline:
        - Galaxia
        - synthetic.perform_pop_syn

    Parameters
    ----------
    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    galactic_config : dict
        Dictionary containing the galactic parameters for the PopSyCLE run

    popsycle_config : dict
        Dictionary containing the PopSyCLE parameters for the PopSyCLE run

    filename_dict : dict
        Dictionary containing the names of the files output by the pipeline

    debugFlag : bool
        If set to True, scripts will be run with a fixed seed that produces
        identical output. If set to False, a random seed will be selected.
        Default is False

    Output
    ------
    None

    """
    # Remove Galaxia output if already exists
    if os.path.exists(filename_dict['ebf_filename']):
        os.remove(filename_dict['ebf_filename'])

    # If debugFlag == True, force Galaxia and PyPopStar (within PopSycLE) to
    # run with fixed seeds
    if debugFlag:
        seed = 0
        set_random_seed = True
    else:
        seed = None
        set_random_seed = False

    # Write out parameters for Galaxia run to disk
    print('-- Generating galaxia params')
    synthetic.write_galaxia_params(
        output_root=output_root,
        longitude=galactic_config['longitude'],
        latitude=galactic_config['latitude'],
        area=galactic_config['area'],
        seed=seed)

    # Run Galaxia from that parameter file
    print('-- Executing galaxia')
    _ = execute('galaxia -r galaxia_params.%s.txt' % output_root)

    # Remove perform_pop_syn output if already exists
    if os.path.exists(filename_dict['hdf5_filename']):
        os.remove(filename_dict['hdf5_filename'])

    # Run perform_pop_syn
    print('-- Executing perform_pop_syn')
    synthetic.perform_pop_syn(
        ebf_file=filename_dict['ebf_filename'],
        output_root=output_root,
        iso_dir=popsycle_config['isochrones_dir'],
        bin_edges_number=popsycle_config['bin_edges_number'],
        BH_kick_speed=popsycle_config['BH_kick_speed'],
        NS_kick_speed=popsycle_config['NS_kick_speed'],
        set_random_seed=set_random_seed)


def run_stage2(output_root,
               popsycle_config, filename_dict,
               parallelFlag=False,
               rank=0,
               comm=None):
    """
    Run stage 2 of the PopSyCLE pipeline:
        - synthetic.calc_events

    Parameters
    ----------
    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    popsycle_config : dict
        Dictionary containing the PopSyCLE parameters for the PopSyCLE run

    filename_dict : dict
        Dictionary containing the names of the files output by the pipeline

    parallelFlag : bool
        If set to True, scripts will be run assuming that it has been launched
        in parallel using mpi4py. If set to False, script will run in serial.
        Default is False

    rank : int
        Number of rank for running script in parallel using mpi4py.
        Default is 0.

    comm : mpi4py.MPI.Intracomm
        Intracommunicator for running script in parallel using mpi4py.
        For example:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD


    Output
    ------
    None

    """
    # Remove calc_events output if already exists
    if os.path.exists(filename_dict['events_filename']):
        os.remove(filename_dict['events_filename'])
    if os.path.exists(filename_dict['blends_filename']):
        os.remove(filename_dict['blends_filename'])

    # Run calc_events
    if rank == 0:
        print('-- Executing calc_events')
    synthetic.calc_events(hdf5_file=filename_dict['hdf5_filename'],
                          output_root2=output_root,
                          radius_cut=popsycle_config['radius_cut'],
                          obs_time=popsycle_config['obs_time'],
                          n_obs=popsycle_config['n_obs'],
                          theta_frac=popsycle_config['theta_frac'],
                          blend_rad=popsycle_config['blend_rad'],
                          overwrite=True)

    # If script is run in parallel, wait for all processes to
    # finish calc_events
    if parallelFlag:
        comm.Barrier()

    if rank == 0:
        # Write a fle to disk stating that there are no events if
        # calc_events does not produce an events file
        if not os.path.exists(filename_dict['events_filename']):
            Path(filename_dict['noevents_filename']).touch()


def run_stage3(output_root, popsycle_config, filename_dict):
    """
    Run stage 3 of the PopSyCLE pipeline:
        - synthetic.refine_events

    Parameters
    ----------
    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    popsycle_config : dict
        Dictionary containing the PopSyCLE parameters for the PopSyCLE run

    filename_dict : dict
        Dictionary containing the names of the files output by the pipeline


    Output
    ------
    None

    """
    # If calc_events failed to produce an event files and instead created a
    # no-events file, exit before running refine_events
    if os.path.exists(filename_dict['noevents_filename']):
        print('No events present, skipping refine_events')
        os.exit(1)

    # Run refine_events
    print('-- Executing refine_events')
    synthetic.refine_events(input_root=output_root,
                            filter_name=popsycle_config['filter_name'],
                            red_law=popsycle_config['red_law'],
                            overwrite=True,
                            output_file='default')


def run():
    description_str = """
    Run the stages of the PopSyCLE pipeline. This executable is intended to be 
    run by slurm scripts generated by `generate_slurm_scripts`. 
    
    Serial stages should be run with 
    `python run_on_slurm.py ...`. Parallel stages should be run with 
    `mpiexec -n 4 run_on_slurm.py ...` or the equivalent way to launch a 
    python script that will run with mpi4py.
    
    Script must be executed in a folder containing a galactic_config file 
    generated by `generate_galactic_config_file`.
    
    Stage 1: (serial)
        - Galaxia
        - synthetic.perform_pop_syn
    Stage 2: (parallel)
        - synthetic.calc_events
    Stage 3: (serial)
        - synthetic.refine_events
        """
    parser = argparse.ArgumentParser(description=description_str)
    parser.add_argument('--output-root', type=str,
                        help='Base filename of the output files',
                        required=True)
    parser.add_argument('--popsycle-config-filename', type=str,
                        help='Name of popsycle_config.yaml file containing '
                             'the PopSyCLE parameters',
                        required=True)
    parser.add_argument('--stage', type=int,
                        help='Number of the stage to be executed. Must be '
                             'either 1, 2 or 3',
                        required=True)
    parser.add_argument('--debug', help='Force Galaxia and PyPopStar '
                                        '(within PopSyCLE) to fix their '
                                        'random seeds to set numbers. This '
                                        'guarantees identical output for '
                                        'all stages.',
                        action='store_true')
    args = parser.parse_args()

    # Load the config files for galactic parameters
    galactic_config = load_galactic_config(args.output_root)
    # Load the config files for popsycle parameters. If the `bin_edges_number`
    # has been set to the string `None`, instead set it to the boolean None.
    popsycle_config = load_popsycle_config(args.popsycle_config_filename)
    if popsycle_config['bin_edges_number'] == 'None':
        popsycle_config['bin_edges_number'] = None

    # Create an isochrones mirror in the current directory
    isochrones_dir = './isochrones'
    if not os.path.exists(isochrones_dir):
        os.symlink(popsycle_config['isochrones_dir'], isochrones_dir)

    # Return the dictionary containing PopSyCLE output filenames
    filename_dict = return_filename_dict(galactic_config['output_root'])

    # Detect parallel processes
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    if rank == 0:
        print('**** Processing stage %i for %s ****' % (args.stage,
                                                        args.output_root))

    # Confirm that stages 1 and 3 are only run with a single processer
    if args.stage in [1, 3] and size != 1:
        print('Stage 1 or 3 must be run with only one rank. Exiting...')
        os.exit(1)

    # Run a different stage depending on the stage parameter
    if args.stage == 1:
        run_stage1(args.output_root,
                   galactic_config, popsycle_config, filename_dict,
                   debugFlag=args.debug)
    elif args.stage == 2:
        run_stage2(args.output_root,
                   popsycle_config, filename_dict,
                   parallelFlag=True, rank=rank, comm=comm)
    elif args.stage == 3:
        run_stage3(args.output_root,
                   popsycle_config, filename_dict)
    else:
        print('Error: stage must be either 1, 2, or 3. Exiting...')
        os.exit(1)


if __name__ == '__main__':
    run()
