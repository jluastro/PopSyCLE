#! /usr/bin/env python
"""
slurm.py
Generate and execute scripts for parallel execution of PopSyCLE runs. Scripts
created will be formatted for submission to a SLURM scheduler.
"""

import subprocess
import yaml
import os


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


def generate_field_config_file(path_run, output_root,
                               longitude, latitude, area):
    """
    Generate a yaml file that contains the microlensing parameters specific to
    a PopSyCLE field.

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
        If set to True, removes all random sampling and forces identical
        output for Galaxia, PyPopStar and PopSyCLE.
        Default False.

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
    # Maximum number of ores per node
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_nodes_max = slurm_config[resource]['n_nodes_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']
    # Get filepath of the run_on_slurm file
    run_filepath = os.path.dirpath(__file__)
    # Create jobname
    jobname = '%s_s%i' % (jobname_base, stage)

    # Template for writing slurm script. Text must be left adjusted.
    mpi_template = """#!/bin/sh
# Job name
#SBATCH --account={account}
#SBATCH --qos={queue}
#SBATCH --constraint={resource}
#SBATCH --nodes={N_nodes}
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
srun -N{N_nodes} -n{N_cores} {path_python} {run_filepath}run.py --output-root={output_root} --stage={stage} --popsycle-config-filename={popsycle_config_filename} {debug_cmd} 
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
    generate_field_config_file(path_run, output_root,
                               longitude, latitude, area)

    # Load the slurm configuration file
    slurm_config = load_slurm_config(slurm_config_filename)

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


def load_slurm_config(slurm_config_filename):
    """
    Load the slurm parameters from the yaml file

    Parameters
    ----------
    slurm_config_filename : str
        Name of slurm_config.yaml file containing the slurm parameters
        that will be used the generate the slurm script header.

    Output
    ------
    slurm_config : dict
        Dictionary containing the slurm parameters for the PopSyCLE run

    """
    # Load the configuration file containing slurm parameters
    with open(slurm_config_filename, 'r') as f:
        slurm_config = yaml.safe_load(f)
    return slurm_config
