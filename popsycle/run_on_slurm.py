#! /usr/bin/env python
"""
run_on_slurm.py
Generate and execute scripts for parallel execution of PopSyCLE runs.
Scripts created will be formatted for submission to a SLURM scheduler.
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


# def submit_script(slurm_config, stage, previous_slurm_job_id=None):
def submit_script(stage, slurm_config, path_run,
                  longitude, latitude, area, walltime,
                  N_nodes_calc_events=None, N_cores_calc_events=None,
                  previous_slurm_job_id=None):
    # Path to the python executable
    path_python = slurm_config['path_python']
    # Project account name to charge
    account = slurm_config['account']
    # Queue
    queue = slurm_config['queue']
    # Number of nodes per run

    # Defult walltime
    if stage == 1:
        n_nodes = 1
        n_cores = 1
    elif stage == 2:
        n_nodes = N_nodes_calc_events
        n_cores = N_cores_calc_events
    elif stage == 3:
        n_nodes = 1
        n_cores = 1
    else:
        raise Exception('stage must be one of [1, 2, 3]')
    # Name of the resource that will be used for the run
    resource = slurm_config['resource']
    # Number of cores per node to use
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_node_max = slurm_config[resource]['n_node_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']

    # Get current directory
    popsycle_directory = os.path.abspath(__file__)

    mpi_template = """#!/bin/sh
    # Job name
    #SBATCH --account={account}
    #SBATCH --qos={queue}
    #SBATCH --constraint={resource}
    #SBATCH --nodes={n_nodes}
    #SBATCH --time={walltime}

    echo "---------------------------"
    date
    echo "Job id = $SLURM_JOBID"
    echo "Proc id = $SLURM_PROCID"
    hostname
    echo "---------------------------"

    cd {path_run}
    srun {dependency} -N{n_nodes} -n{n_cores} {path_python} {popsycle_directory}/run.py --stage={stage} 

    date
    echo
    "All done!"
    """

    # Check that the specified number of nodes does not exceed the resource max
    if n_nodes > n_node_max:
        print('Error: specified number of nodes exceeds limit. Exiting')
        os.exit()

    # Make a run directory within the field file
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # Add a dependency if previous_slurm_id is not none
    if previous_slurm_job_id is None:
        dependency = ''
    else:
        dependency = '--dependency=afterok:{0}'.format(previous_slurm_job_id)

    # Get the total number of cores
    n_cores = n_nodes * n_cores_per_node

    # Populate the mpi_template specified inputs
    job_script = mpi_template.format(**locals())

    # Save the script
    if not path_run.endswith('/'):
        path_run = path_run + '/'
    script_filename = path_run + 'run_popsycle_{0}.sh'.format(stage)
    with open(script_filename, 'w') as f:
        f.write(job_script)

    # Submit the job
    stdout, stderr = execute('sbatch {0}'.format(script_filename))

    print('Submitted job {0} to {1}'.format(script_filename, resource))

    previous_slurm_job_id = stdout.decode().replace('\n', '').split()[-1]

    return previous_slurm_job_id


def submit_pipeline(slurm_config_file, path_run, longitude, latitude, area,
                    N_nodes_calc_events, N_cores_calc_events,
                    walltime_stage1, walltime_stage2,
                    walltime_stage3):
    with open(slurm_config_file, 'r') as stream:
        slurm_config = yaml.safe_load(stream)

    slurm_job_id1 = submit_script(1, slurm_config, path_run,
                                  longitude, latitude, area, walltime_stage1)
    slurm_job_id2 = submit_script(2, slurm_config, path_run, slurm_job_id1,
                                  longitude, latitude, area, walltime_stage2,
                                  N_nodes_calc_events=N_nodes_calc_events,
                                  N_cores_calc_events=N_cores_calc_events,
                                  previous_slurm_job_id=slurm_job_id1)
    _ = submit_script(3, slurm_config, path_run, slurm_job_id2,
                      longitude, latitude, area, walltime_stage3,
                      previous_slurm_job_id=slurm_job_id2)
