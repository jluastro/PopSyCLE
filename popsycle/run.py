#! /usr/bin/env python
"""
run.py
Executable to run the PopSyCLE pipeline in stages.
"""

import yaml
import os
from pathlib import Path
import argparse
from argparse import RawTextHelpFormatter
from popsycle import synthetic
from popsycle.slurm import execute


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
               filename_dict, debug=False):
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

    debug : bool
        If set to True, removes all random sampling and forces identical
        output for Galaxia, PyPopStar and PopSyCLE.
        Default False.

    Output
    ------
    None

    """
    # Remove Galaxia output if already exists
    if os.path.exists(filename_dict['ebf_filename']):
        os.remove(filename_dict['ebf_filename'])

    # If debugFlag == True, force Galaxia and PyPopStar (within PopSycLE) to
    # run with fixed seeds
    if debug:
        seed = 0
    else:
        seed = None

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
        debug=debug)


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
        in parallel using mpi4py. If set to False, script will run in serial
        and skip over mpi4py import.
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

    # Run calc_events
    if rank == 0:
        print('-- Executing calc_events')
        # Remove calc_events output if already exists
        if os.path.exists(filename_dict['events_filename']):
            os.remove(filename_dict['events_filename'])
        if os.path.exists(filename_dict['blends_filename']):
            os.remove(filename_dict['blends_filename'])

    # If execution is parallel, all processes wait to proceed
    if parallelFlag:
        comm.Barrier()

    synthetic.calc_events(hdf5_file=filename_dict['hdf5_filename'],
                          output_root2=output_root,
                          radius_cut=popsycle_config['radius_cut'],
                          obs_time=popsycle_config['obs_time'],
                          n_obs=popsycle_config['n_obs'],
                          theta_frac=popsycle_config['theta_frac'],
                          blend_rad=popsycle_config['blend_rad'],
                          parallelFlag=True,
                          overwrite=True)

    # If execution is parallel, wait for all processes to finish calc_events
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
    Run the stages of the PopSyCLE pipeline. This executable can be either 
    run by slurm scripts generated by `generate_slurm_scripts` or from the 
    command line.

    Serial stages should be run with 
    `python {PATH_TO_FILE}/run.py ...`. Parallel stages should be run with 
    `mpiexec -n 4 {PATH_TO_FILE}/run.py ...` or the equivalent way to launch a 
    python script that will run with mpi4py.

    Script must be executed in a folder containing a galactic_config file 
    generated by `popsycle.slurm.generate_galactic_config_file`.

    Stage 1: Serial
        - Galaxia
        - synthetic.perform_pop_syn
    Stage 2: Parallel
        - synthetic.calc_events
    Stage 3: Serial
        - synthetic.refine_events
        """

    parser = argparse.ArgumentParser(description=description_str,
                                     formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group(title='Required')
    required.add_argument('--stage', type=int,
                        help='Number of the stage to be executed. Must be '
                             'either 1, 2 or 3.',
                        required=True)
    required.add_argument('--output-root', type=str,
                        help='Base filename of the output files. '
                             'Default: root0',
                        default='root0')
    required.add_argument('--popsycle-config-filename', type=str,
                        help='Name of popsycle_config.yaml file containing '
                             'the PopSyCLE parameters. '
                             'Default: popsycle_config.yaml',
                        default='popsycle_config.yaml')

    optional = parser.add_argument_group(title='Optional')
    optional.add_argument('--debug', help='Force Galaxia and PyPopStar '
                                        '(within PopSyCLE) to fix their '
                                        'random seeds to set numbers. This '
                                        'guarantees identical output for '
                                        'all stages.',
                        action='store_true')
    optional.add_argument('--serial', help="Force all stages of the pipeline to "
                                         "run in serial regardless of the "
                                         "presence of multiple cores or how "
                                         "the script is launched. Running "
                                         "run.py with --serial will avoid "
                                         "any 'import mpi4py' statements.",
                        action='store_true')
    args = parser.parse_args()
    print(args)
    import sys
    sys.exit(0)

    # Check for galactic config file. Exit if not present.
    config_filename = 'galactic_config.{0}.yaml'.format(args.output_root)
    if not os.path.exists(config_filename):
        print("""Error: Galactic configuration file {0} missing from current 
        folder, cannot continue. In order to execute run.py, generate a 
        galactic configuration file using 
        popsycle.slurm.generate_galactic_config_file. Exiting...""")
        os.exit(1)

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
    filename_dict = return_filename_dict(args.output_root)

    # Detect parallel processes
    if args.serial:
        comm = None
        rank = 0
        size = 1
        parallelFlag = False
    else:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        parallelFlag = True

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
                   debug=args.debug)
    elif args.stage == 2:
        run_stage2(args.output_root,
                   popsycle_config, filename_dict,
                   parallelFlag=parallelFlag, rank=rank, comm=comm)
    elif args.stage == 3:
        run_stage3(args.output_root,
                   popsycle_config, filename_dict)
    else:
        print('Error: stage must be either 1, 2, or 3. Exiting...')
        os.exit(1)


if __name__ == '__main__':
    run()
