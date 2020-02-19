import numpy as np
import pylab as plt
import pytest
from popsycle import utils

def test_add_precision64_single():
    res1 = utils.add_precision64(1.0, -4)

    assert np.abs(res1 - 1.0) > 1e-6
    
    return

def test_add_precision64_array():
    in_vals = np.ones(10)
    out_vals = utils.add_precision64(in_vals, -4)

    np.testing.assert_almost_equal(out_vals, 1.0, 1)
    
    return

def test_sample_spherical_single_speed_single_samples_3dim():
    # Testing single speed with single sample in default 3 dimensions
    N_samples = 1
    speed_input = 1
    output = utils.sample_spherical(N_samples, speed_input)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == 3
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)


def test_sample_spherical_single_speed_multiple_samples_3dim():
    # Testing single speed with multiple samples in default 3 dimensions
    N_samples = 5
    speed_input = 1
    output = utils.sample_spherical(N_samples, speed_input)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == 3
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)


def test_sample_spherical_multiple_speeds_multiple_samples_3dim():
    # Testing multiple speeds with multiple samples in default 3 dimensions
    N_samples = 5
    speed_input = np.arange(1, 6).astype(float)
    output = utils.sample_spherical(N_samples, speed_input)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == 3
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)


def test_sample_spherical_multiple_speeds_single_sample_3dim():
    # Testing sample_spherical breaks correctly if speed of incorrect
    # length is inserted
    with pytest.raises(ValueError):
        N_samples = 1
        speed_input = np.arange(1, 6).astype(float)
        _ = utils.sample_spherical(N_samples, speed_input)


def test_sample_spherical_single_speed_single_samples_4dim():
    # Testing single speed with single sample in 4 dimensions
    N_samples = 1
    speed_input = 1
    ndim = 4
    output = utils.sample_spherical(N_samples, speed_input, ndim=ndim)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == ndim
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)


def test_sample_spherical_single_speed_multiple_samples_4dim():
    # Testing single speed with multiple samples in 4 dimensions
    N_samples = 5
    speed_input = 1
    ndim = 4
    output = utils.sample_spherical(N_samples, speed_input, ndim=ndim)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == ndim
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)


def test_sample_spherical_multiple_speeds_multiple_samples_4dim():
    # Testing multiple speeds with multiple samples in 4 dimensions
    N_samples = 5
    speed_input = np.arange(1, 6).astype(float)
    ndim = 4
    output = utils.sample_spherical(N_samples, speed_input, ndim=ndim)
    speed_output = np.sqrt(np.sum(output ** 2, axis=0))

    assert output.shape[0] == ndim
    assert output.shape[1] == N_samples
    np.testing.assert_almost_equal(speed_input, speed_output, decimal=7)