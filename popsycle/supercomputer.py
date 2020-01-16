#! /usr/bin/env python
"""
supercomputer.py
Generate and execute slurm scripts for parallel execution of PopSyCLE runs.
"""

import subprocess

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