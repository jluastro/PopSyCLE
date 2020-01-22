#! /usr/bin/env python
"""
slurm_launcher.py
Generate and execute scripts for parallel execution of PopSyCLE runs.
Scripts created will be formatted for submission to a SLURM scheduler.
"""

import subprocess
from argparse import ArgumentParser
import yaml
import os


def write_galaxia_params(output_dir, output_basename='ZTF1',
                         seed=0,
                         longitude=45.19260648,
                         latitude=4.93717557,
                         area=10.00):
    params = [
        "outputFile %s" % output_basename,
        "outputDir %s" % output_dir,
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
        "surveyArea %.2f" % area,
        "fSample 1",
        "popID -1",
        "warpFlareOn 1",
        "seed %i" % seed,
        "r_max 30",
        "starType 0",
        "photoError 0"
    ]
    with open(output_dir + '/galaxia_params.%i.txt' % seed, 'w') as f:
        for param in params:
            f.write(param + '\n')


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

    # Get current directory
    popsycle_directory = os.path.abspath(__file__)

    mpi_template = """#!/bin/sh
    # Job name
    #SBATCH --account={account}
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
    if previous_slurm_id:
        dependency = '--dependency=afterok:{previous_slurm_id}'
    else:
        dependency = ''

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

    print
    'Submitted job {0} to {1}'.format(script_file, resource)

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
