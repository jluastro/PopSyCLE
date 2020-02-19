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
