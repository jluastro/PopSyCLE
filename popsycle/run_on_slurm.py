#! /usr/bin/env python
"""
run_on_slurm.py
Generate and execute scripts for parallel execution of PopSyCLE runs.
Scripts created will be formatted for submission to a SLURM scheduler.
"""

import subprocess
import yaml
import os
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


def generate_microlensing_param_file(path_run, output_root,
                                     longitude, latitude, area):
    param_filename = '{0}/microlensing_params.{1}.yaml'.format(path_run,
                                                               output_root)
    with open(param_filename, 'w') as f:
        f.write('longitude: %f\n' % longitude)
        f.write('latitude: %f\n' % latitude)
        f.write('area: %f\n' % area)
        f.write('output_root: %s\n' % output_root)


# def submit_script(slurm_config, stage, previous_slurm_job_id=None):
def submit_script(stage, slurm_config, microlensing_config_filename,
                  path_run, output_root,
                  jobname_base, walltime,
                  N_nodes_calc_events=None, N_cores_calc_events=None,
                  previous_slurm_job_id=None,
                  skipSubmitFlag=False, debugFlag=False):
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
        raise Exception('Error: stage must be one of [1, 2, 3]. Exiting...')
    # Name of the resource that will be used for the run
    resource = slurm_config['resource']
    # Number of cores per node to use
    n_cores_per_node = slurm_config[resource]['n_cores_per_node']
    # Maximum number of nodes
    n_nodes_max = slurm_config[resource]['n_nodes_max']
    # Maximum walltime (hours)
    walltime_max = slurm_config[resource]['walltime_max']
    # Get current directory
    run_on_slurm_path = os.path.abspath(__file__)
    # Create jobname
    jobname = '%s_s%i' % (jobname_base, stage)

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
cd {path_run}
srun -N{n_nodes} -n{n_cores} {path_python} {run_on_slurm_path}π --output-root={output_root} --stage={stage} --microlensing-config-filename={microlensing_config_filename} 
date
echo
"All done!"
"""

    # Check that the specified number of nodes does not exceed the resource max
    if n_nodes > n_nodes_max:
        raise Exception('Error: specified number of nodes exceeds limit. '
                        'Exiting...')
    if n_cores > n_cores_per_node:
        raise Exception('Error: specified number of cores exceeds limit. '
                        'Exiting...')

    # Make a run directory within the field file
    if not os.path.exists(path_run):
        os.makedirs(path_run)

    # Populate the mpi_template specified inputs
    job_script = mpi_template.format(**locals())

    # Save the script
    if not path_run.endswith('/'):
        path_run = path_run + '/'
    script_filename = path_run + 'run_popsycle_{0}.sh'.format(stage)
    with open(script_filename, 'w') as f:
        f.write(job_script)

    # Submit the job
    # Add a dependency if previous_slurm_id is not none
    if skipSubmitFlag:
        previous_slurm_job_id = None
    else:
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

        previous_slurm_job_id = stdout.decode().replace('\n', '').split()[-1]

    return previous_slurm_job_id


def submit_pipeline(slurm_config_file, microlensing_config_filename,
                    path_run, output_root,
                    longitude, latitude, area,
                    N_nodes_calc_events, N_cores_calc_events,
                    walltime_stage1, walltime_stage2,
                    walltime_stage3, skipSubmitFlag=False, debugFlag=False):
    generate_microlensing_param_file(path_run, output_root,
                                     longitude, latitude, area)

    with open(slurm_config_file, 'r') as f:
        slurm_config = yaml.safe_load(f)

    jobname_base = 'l%.1f_b%.1f' % (longitude, latitude)
    slurm_job_id1 = submit_script(1, slurm_config,
                                  microlensing_config_filename,
                                  path_run, output_root,
                                  jobname_base, walltime_stage1,
                                  skipSubmitFlag=skipSubmitFlag,
                                  debugFlag=debugFlag)
    slurm_job_id2 = submit_script(2, slurm_config,
                                  microlensing_config_filename,
                                  path_run, output_root,
                                  jobname_base, walltime_stage2,
                                  N_nodes_calc_events=N_nodes_calc_events,
                                  N_cores_calc_events=N_cores_calc_events,
                                  previous_slurm_job_id=slurm_job_id1,
                                  skipSubmitFlag=skipSubmitFlag,
                                  debugFlag=debugFlag)
    _ = submit_script(3, slurm_config, microlensing_config_filename,
                      path_run, output_root,
                      jobname_base, walltime_stage3,
                      previous_slurm_job_id=slurm_job_id2,
                      skipSubmitFlag=skipSubmitFlag,
                      debugFlag=debugFlag)


def load_microlensing_params(output_root):
    params_filename = 'microlensing_params.{0}.yaml'.format(output_root)
    with open(params_filename, 'r') as f:
        microlensing_params = yaml.safe_load(f)
    return microlensing_params


def load_microlensing_config(microlensing_config_filename):
    with open(microlensing_config_filename, 'r') as f:
        microlensing_config = yaml.safe_load(f)
    return microlensing_config


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-root', type=str,
                        required=True)
    parser.add_argument('--microlensing-config-filename', type=str,
                        required=True)
    parser.add_argument('--stage', type=int,
                        required=True)
    args = parser.parse_args()

    microlensing_params = load_microlensing_params(args.output_root)
    microlensing_config = load_microlensing_config(
        args.microlensing_config_filename)
    if args.stage == 1:
        # Galaxia
        random_seed = synthetic.write_galaxia_params(
            output_root=microlensing_params['output_root'],
            longitude=microlensing_params['longitude'],
            latitude=microlensing_params['latitude'],
            area=microlensing_params['area'],
            seed=None)
        _ = execute('galaxia -r galaxia_params.%i.txt' % random_seed)
        # perform_pop_syn
        ebf_filename = '%s.ebf' % microlensing_params['output_root']
        isochrones_dir = './isochrones'
        if not os.path.exists(isochrones_dir):
            os.symlink(microlensing_config['isochrones_dir'], isochrones_dir)

        synthetic.perform_pop_syn(
            ebf_file=ebf_filename,
            output_root=microlensing_params['output_root'],
            iso_dir=microlensing_config['isochrones_dir'],
            bin_edges_number=microlensing_config['bin_edges_number'],
            BH_kick_speed=microlensing_config['BH_kick_speed'],
            NS_kick_speed=microlensing_config['NS_kick_speed'],
            overwrite=True)
    elif args.stage == 2:
        # calc_events
        hdf5_file = '%s.h5' % microlensing_params['output_root']
        synthetic.calc_events(hdf5_file=hdf5_file,
                              output_root2=microlensing_params['output_root'],
                              radius_cut=2,
                              obs_time=1095,
                              n_obs=5,
                              theta_frac=2,
                              blend_rad=0.75,
                              overwrite=True)
    elif args.stage == 3:
        synthetic.refine_events(input_root=microlensing_params['output_root'],
                                filter_name=microlensing_config['filter_name'],
                                red_law=microlensing_config['red_law'],
                                overwrite=True,
                                output_file='default')
    else:
        raise Exception('Error: stage must be one of [1, 2, 3]. Exiting...')


if __name__ == '__main__':
    run()
