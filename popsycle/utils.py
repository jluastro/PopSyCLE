import numpy as np
import h5py
import math
from astropy import units
from scipy.stats import maxwell
import astropy.coordinates as coord
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.table import Table
from astropy.table import vstack
from popstar.imf import imf
from popstar import synthetic, evolution, reddening, ifmr
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
import time
import datetime
from popsycle import ebf
import gc
import subprocess
import os
from sklearn import neighbors
import itertools
from multiprocessing import Pool
import yaml
import inspect
from popstar import atmospheres
from popstar.imf import multiplicity
from scipy.interpolate import griddata
import numpy.lib.recfunctions as rfn


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


