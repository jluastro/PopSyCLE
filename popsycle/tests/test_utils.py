import numpy as np
import pylab as plt
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
