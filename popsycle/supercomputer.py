#! /usr/bin/env python
"""
supercomputer.py
Generate and execute slurm scripts for parallel execution of PopSyCLE runs.
Scripts created will be formatted for submission to a SLURM scheduler.
"""

import subprocess
from argparse import ArgumentParser
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


def generate_script(args,
                    stage,
                    previous_slurm_id=None):
    slurm_config = yaml.load(open(args.slurm_config_file))

    # Folder for outputting PopSyCLE data
    path_run = slurm_config['path_run']
    # Path to the python executable
    path_python = slurm_config['path_python']
    # Path to the pipeline
    path_pipeline = slurm_config['path_pipeline']
    # Project account name to charge
    account = slurm_config['account']
    # Flag to mute warning output to log file
    mute_warnings = slurm_config['mute_warnings']
    # MPI version to use
    version_mpi = slurm_config['version_mpi']
    # HDF5 version to use
    version_hdf5 = slurm_config['version_hdf5']
    # Queue
    queue = slurm_config['queue']
    # Number of nodes per run 
    n_nodes = slurm_config['n_nodes']
    # Defult walltime
    if stage == 1:
        walltime = args.walltime_1
    elif stage == 2:
        walltime = args.walltime_2
    elif stage == 3:
        walltime = args.walltime_3
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

    mpi_template = """#!/bin/csh
    # Job name
    #SBATCH -N {name}
    #SBATCH -A {account}
    #
    # Combine stdout and stderr
    #SBATCH -j oe
    #SBATCH -o {path_run}/slurm.log
    #SBATCH -m bea
    #
    #SBATCH -q {queue}
    #SBATCH -l nodes={n_nodes}
    #SBATCH -l walltime={walltime}
    #SBATCH -V

    echo "---------------------------"
    date
    echo "Job id = $SLURM_JOBID"
    echo "Proc id = $SLURM_PROCID"
    hostname

    # Specify a consistent mpi version
    use {version_mpi}
    # Specify a consistent hdf5 version
    use {version_hdf5}

    # There is no need to activate the environment as long as you specify the
    # environments python.
    """

    if previous_slurm_id:
        mpi_template_append = """
        srun --dependency=afterok:{previous_slurm_id} -N{n_nodes} -n{n_cores} {path_python} {warning_opt} {path_pipeline} -f {file_name} -c {inject_params}

        date
        echo
        "All done!"
        """
    else:
        mpi_template_append = """
        srun -N{n_nodes} -n{n_cores} {path_python} {warning_opt} {path_pipeline} -f {file_name} -c {inject_params}

        date
        echo
        "All done!"
        """
    mpi_template += mpi_template_append

    # Check that the specified number of nodes does not exceed the resource max
    if n_nodes > n_node_max:
        print 'Error: specified number of nodes exceeds limit. Exiting'
        os.exit()

    # Make a run directory within the field file
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # Set the warning log option for output log files
    if mute_warnings == True:
        warning_opt = '-W ignore'
    else:
        warning_opt = ''

    # Get the total number of cores
    n_cores = n_nodes * n_cores_per_node

    # Populate the mpi_template specified inputs
    job_script = mpi_template.format(**locals())

    # Save the script
    if not path_run.endswith('/'):
        path_run = path_run + '/'
    script_file = path_run + 'run_popsycle_{0}.sh'.format(stage)
    with open(script_file, 'w') as f:
        f.write(job_script)

    # Submit the job
    stdout, stderr = execute('sbatch {0}'.format(script_file))

    print 'Submitted job {0} to {1}'.format(script_file, resource)

    slurm_job_id = stdout[0]  # FIXME Implement correct formatting

    return slurm_job_id


def main():
    parser = ArgumentParser(description='Generate and submit slurm scripts '
                                        'for running PopSyCle')
    arguments = parser.add_argument_group('arguments')
    arguments.add_argument('--galaxy-l', type=float,
                           help='Galactic longitude (degrees)',
                           required=True)
    arguments.add_argument('--galaxy-b', type=float,
                           help='Galactic latitude (degrees)',
                           required=True)
    arguments.add_argument('--area', type=float,
                           help='Area on sky (square degrees)',
                           required=True)
    arguments.add_argument('--N-nodes', type=int,
                           help='Number of nodes for running calc_events',
                           required=True)
    arguments.add_argument('--walltime-1', type=float,
                           help='Walltime (hours) for running Galaxia '
                                'and perform_pop_sny',
                           required=True)
    arguments.add_argument('--walltime-2', type=float,
                           help='Walltime (hours) for running calc_events',
                           required=True)
    arguments.add_argument('--walltime-3', type=float,
                           help='Walltime (hours) for running refine_events'
                                'and perform_pop_sny',
                           required=True)
    arguments.add_argument('--slurm-config-file', type=str,
                           default='slurm-config.yaml',
                           help='A user created YAML configuration file that '
                                'specifies the scheduler inputs')

    args = parser.parse_args()

    slurm_job_id1 = generate_script(args=args,
                                    stage=1)
    slurm_job_id2 = generate_script(args=args,
                                    stage=2,
                                    previous_slurm_id=slurm_job_id1)
    slurm_job_id3 = generate_script(args=args,
                                    stage=3,
                                    previous_slurm_id=slurm_job_id2)


if __name__ == '__main__':
    main()
