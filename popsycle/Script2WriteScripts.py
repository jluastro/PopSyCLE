#!/usr/bin/env python
"""
This code will run the likelihood detection pipeline on all the tile files in
the specified field folder. It is designed to be run on a Livermore Computing
Moab and Slurm system (e.g. quartz or cab).
"""
from __future__ import division
from argparse import ArgumentParser
import yaml
import os
import glob
import datetime
import numpy as np

mpi_template = """#!/bin/csh
# Job name
#MSUB -N {name}
#MSUB -A {account}
#
# Combine stdout and stderr
#MSUB -j oe
#MSUB -o {path_run}logs/{name}_detect.log
#MSUB -m bea
#
#MSUB -q {queue}
#MSUB -l nodes={n_nodes}
#MSUB -l walltime={walltime}
#MSUB -V

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
# environements python.
srun -N{n_nodes} -n{n_cores} {path_python} {warning_opt} {path_pipeline} -f {file_name} -c {inject_params}

date
echo "All done!"
"""

def get_walltime(n_cores, file_size):
    '''

    # Convert file size to MB for scaling relation input

    '''
    file_size = file_size / 1e6

    minutes = np.int(np.round(file_size*1.2+60))

    return minutes


if __name__ == '__main__':
    # Inatilize the parser with the default config file.
    parser = ArgumentParser(description=
                            ('A program to submit a moab script for each tile ',
                             'file in the specified field directory.'))
    parser.add_argument('-c', '--configfile', type=str,
                        default='run_field_config.yaml',
                        help=('A user created YAML configuration file that ',
                              'specifies the parameter inputs.'))
    parser.add_argument('-f', '--path_data', type=str,
                        default=None,
                        help=('Optional command line argument to specify the ',
                              'field directory path. By defult this will use ',
                              'the path specified in the configfile.'))
    # read in the arguments
    args = parser.parse_args()
    # Read the yaml configuraton file.
    params = yaml.load(open(args.configfile))

    # Assign the user specified parameters
    # Overwrite the config file param if command line argument provided
    if args.path_data == None:
        path_data = params['path_data']
    else:
        path_data = args.path_data
    # Path to the python executable
    path_python = params['path_python']
    # Path to the pipeline
    path_pipeline = params['path_pipeline']
    # Project account name to charge
    account = params['account']
    # Flag to mute warning output to log file
    mute_warnings = params['mute_warnings']
    # MPI version to use
    version_mpi = params['version_mpi']
    # HDF5 version to use
    version_hdf5 = params['version_hdf5']
    # Queue
    queue = params['queue']
    # Number of nodes per run 
    n_nodes = params['n_nodes']
    # Defult walltime
    walltime_default = params['walltime_default']
    # Name of the resource that will be used for the run
    resource = params['resource']
    # Number of cores per node to use
    n_cores_per_node = params[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_node_max = params[resource]['n_node_max']
    # Maximum walltime (hours)
    walltime_max = params[resource]['walltime_max']

    # Check that the specified number of nodes does not exceed the resource max
    if n_nodes > n_node_max:
        print 'Error: specified number of nodes exceeds limit. Exiting'
        os.exit()

    # Make a run directory within the field file
    path_run = path_data + '/micro_detect/'
    if not os.path.exists(path_run):
        os.makedirs(path_run)
    # Make a log directory within the run directory
    path_log = path_run + '/logs/'
    if not os.path.exists(path_log):
        os.makedirs(path_log)

    # Set the warning log option for output log files
    if mute_warnings == True:
        warning_opt = '-W ignore'
    else:
        warning_opt = ''

    # Get the total number of cores
    n_cores = n_nodes * n_cores_per_node

    #specifcy the inject parameter yaml
    inject_params = 'params_MW.yaml'

    # Determine the list of tile files in the directory
    file_list = glob.glob(path_data+'/*.gz')
    # Loop through the file list submitting a seperate job for each file
    for index, file_name in enumerate(file_list):
        # Get the file size
        file_size = os.path.getsize(file_name)
        if file_size < 1000:
            # The the file is likely empty or contains very few stars.
            print 'Warning: skipping {0}. File size = {1} kB and is likely empty or contains very few stars.'.format(file_name, file_size / 1000.0)
            continue
        # get the basename of the file without path or extension
        name = os.path.basename(os.path.splitext(file_name)[0])

        # Estimate the walltime given the file size
        if walltime_default == 'auto':
            if resource != 'quartz':
                print 'Warning: empirical estimates for walltime scaleing currently only implemented for quartz. Walltime estimates may be highly inaccurate.'
            t_minutes = get_walltime(n_cores, file_size)
            # Verify that the estimated walltime does not exceed the limits
            if t_minutes / 60.0 > walltime_max:
                print 'Warning: The estimated walltime of {0:0.1f} hours exceeds the limit of {1}. Will skip tile {2}.'.format(t_minutes / 60.0, walltime_max, file_name)
                continue
            # Convert the estimated walltime to the correct format
            t = datetime.timedelta(minutes=t_minutes)
            walltime = str(t)
        else:
            # Use the same user specified walltime for all runs
            walltime = walltime_default

        # Populate the mpi_template specified inputs
        job_script = mpi_template.format(**locals())

        # Save the script 
        script_file = path_run + '{0}.sh'.format(name)
        f = open(script_file, 'w')
        f.write(job_script)
        f.close()

        # Submit the job
        os.system('msub {0}'.format(script_file))

        print 'Submitted job {0} to {1}'.format(name, resource)
