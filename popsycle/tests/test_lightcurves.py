import pytest
import os
from popsycle import lightcurves
from astropy.table import Table
from bagle import model


def test_get_bagle_model_list():
    test_filepath = os.path.dirname(__file__)

    # Specify photometric system variables.
    photom_sys = 'ubv'
    filter_name = 'I'
    red_law = 'Damineli16'
    n_multi_proc = 6

    in_root = f'data_test/test_Mrun_refined_events_{photom_sys}_{filter_name}_{red_law}'

    full_in_root = test_filepath + '/' + in_root

    events_tab = Table.read(full_in_root + '_rb.fits')
    comps_tab = Table.read(full_in_root + '_companions_rb.fits')
    lcurves_tab = Table.read(full_in_root + '_rb_lightcurves.fits')

    model_list = lightcurves.get_bagle_model_list(events_tab, comps_tab, lcurves_tab,
                                                  photom_sys, filter_name, red_law,
                                                  n_multi_proc=n_multi_proc)

    assert isinstance(model_list[0], model.PSBL_PhotAstrom_Par_EllOrbs_Param7)
    assert len(model_list) == len(events_tab)

    return


#FIXME Add test for lightcurve production
def test_psbl_multi_lightcurve():
    pass

def test_

