import pytest
import os
from popsycle import analysis

def test_get_star_system_pos_mag():
    test_filepath = os.path.dirname(__file__)
    in_root = 'data_test/test_correct_Mrun'

    full_in_root = test_filepath + '/' + in_root

    h5_file = full_in_root + '.h5'

    df = analysis.get_star_system_pos_mag(h5_file)

    assert 'm_ubv_I_app' in df.columns
    assert 'glat' in df.columns
    assert 'x' not in df.columns

    return

def test_count_stars_hdf5():
    test_filepath = os.path.dirname(__file__)
    in_root = 'data_test/test_correct_Mrun'

    full_in_root = test_filepath + '/' + in_root

    h5_file = full_in_root + '.h5'

    N_stars = analysis.count_stars_hdf5(h5_file, mag_threshold=100)
    assert N_stars == 37327

    N_stars = analysis.count_stars_hdf5(h5_file, mag_threshold=21)
    assert N_stars == 890

    return

