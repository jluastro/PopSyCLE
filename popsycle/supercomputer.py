#! /usr/bin/env python
"""
supercomputer.py
Generate and execute slurm scripts for parallel execution of PopSyCLE runs.
"""

import subprocess
import argparse


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


def main():
    parser = argparse.ArgumentParser()
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
    arguments.add_argument('--N-cores-per-node', type=int,
                           help='Number of cores on a single node',
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

    args = parser.parse_args()


if __name__ == '__main__':
    main()
