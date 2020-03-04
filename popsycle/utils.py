#! /usr/bin/env python
"""
utils.py
Functions (and their associated functions) for utilities.
"""
import subprocess
import numpy as np


def add_precision64(input_array, power):
    """
    Need more precision for kdtree to run properly. Convert inputs from
    float32 to float64, and add a random perturbation beginning in the
    nths decimal place.

    Parameters
    ----------
    input_array : float or array (float32)
        Thing that needs more precision.

    power : float
        To what place you want the perturbation.

    Return
    ------
    output_array : float or array (float64)
        Thing that has more precision.

    """
    input_array = np.atleast_1d(input_array)

    # Perturb.
    pert = 10 ** power * (np.random.rand(len(input_array)) - 0.5)

    # Convert to float64.
    output_array = np.atleast_1d(np.float64(input_array))

    # Add the perturbation.
    output_array = output_array + pert

    return output_array


def sample_spherical(npoints, speed, ndim=3):
    """
    Randomly sample points on a sphere.
    I found this code on stackexchange.

    Parameters
    ---------
    npoints : float
        The number of points you want to generate.

    speed : float
        The radius of the sphere (aka the magnitude of the vectors.)

    dim : float
        The dimension of the space in which the sphere is embedded
        (ndim = 3 samples points on a 2-sphere, aka a "normal" sphere)

    Return
    ------
    An array of the vectors.
    """
    # Check that the speed vector is either a float or an array of length npoints
    if type(speed) != int:
        if type(speed) != float:
            if len(speed) != npoints:
                raise ValueError("{speed} must be either an int, float "
                                 "or array of length {npoints}")

    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    vec *= speed
    return vec


def return_nearest_gridpoint(grid, x_grid_arr, y_grid_arr, x_data, y_data):
    """
    Algorithm for finding the nearest grid cell on a 2D array given a
    datapoint that falls within the bounds of the 2D array.

    Parameters
    ----------
    grid : 2D numpy array
        2D array with size (len(y_grid_arr), len(x_grid_array))

    x_grid_arr : numpy array
        2D grid indices in the x-dimension

    y_grid_arr : numpy array
        2D grid indices in the y-dimension

    x_data : numpy array
        x-coordinate for data that will be located onto the grid

    y_data : numpy array
        y-coordinate for data that will be located onto the grid

    Output
    ------
    gridpoint_arr : numpy array
        list of nearest cell values on the grid at
        the location of (x_data, y_data)

    """
    # Convert x_data and y_data to array if single data point is received
    x_data = np.atleast_1d(x_data)
    y_data = np.atleast_1d(y_data)

    gridpoint_arr = []
    for x, y in zip(x_data, y_data):
        # Loop through x_data and y_data
        if np.isnan(x) or np.isnan(y):
            # If either x_data or y_data is nan, return nan
            gridpoint = np.nan
        else:
            # Find location on the grid where x_data and y_data are
            # closest to the grid indices
            x_idx = np.argmin(np.abs(x - x_grid_arr))
            y_idx = np.argmin(np.abs(y - y_grid_arr))
            gridpoint = grid[y_idx, x_idx]
        gridpoint_arr.append(gridpoint)

    # Convert gridpoint_arr into numpy array
    gridpoint_arr = np.array(gridpoint_arr)

    # If only a single data point was received, return a single value
    if len(gridpoint_arr) == 1:
        gridpoint_arr = gridpoint_arr[0]

    return gridpoint_arr


def execute(cmd, shell=False):
    """
    Executes a command line instruction, captures the stdout and stderr

    Parameters
    ----------
    cmd : str
        Command line instruction, including any executables and parameters

    Optional Parameters
    -------------------
    shell : bool
        Determines if the command is run through the shell. Default is False.

    Outputs
    -------
    stdout : str
        Contains the standard output of the executed process

    stderr : str
        Contains the standard error of the executed process

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