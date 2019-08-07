from popsycle import synthetic
from astropy.coordinates import SkyCoord # High-level coordinates
import astropy.coordinates as coord
from astropy.table import Table
import ebf
from astropy import units
import numpy as np
import h5py
import pdb

masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcday = 1.0 * (3.086 * 10**16)**-1 * 86400.0
au_to_kpc = 4.848 * 10**-9

def test_calc_blends():
    test_data_dir = '/u/casey/scratch/work/microlens/galaxia_test/OGLE672/'

    hdf5_file = test_data_dir + 'OGLE672.h5'
    event_file = test_data_dir + 'OGLE672_events.fits'

    # Identify the l, b patches we will process (similar to calc_events)
    ll = 10
    bb = 10
    print('Working on loop ll, bb = ', ll, bb)
    name00 = 'l' + str(ll) + 'b' + str(bb)
    name01 = 'l' + str(ll) + 'b' + str(bb + 1)
    name10 = 'l' + str(ll + 1) + 'b' + str(bb)
    name11 = 'l' + str(ll + 1) + 'b' + str(bb + 1)
    print(name00, name01, name10, name11)
    
    hf = h5py.File(hdf5_file, 'r')
    bigpatch = np.concatenate((hf[name00], hf[name01], hf[name10], hf[name11]), axis=1)    
    hf.close()
    print(bigpatch.shape)

    time = 0.0
    r_t = bigpatch[9] + time * bigpatch[12] * kms_to_kpcday # kpc
    b_t = bigpatch[10] + time * bigpatch[13] * masyr_to_degday # deg
    l_t = bigpatch[11] + time * (bigpatch[14] / np.cos(np.radians(bigpatch[10]))) * masyr_to_degday # deg
                
    ##########
    # Define coordinates.
    ##########
    c = SkyCoord(frame = 'galactic', l = l_t * units.deg, b = b_t * units.deg)
    foo = coord.match_coordinates_sky(c, c, nthneighbor = 2)

    # Load up the events table
    event_lbt_tmp = Table.read(event_file)

    # Trim it down for faster processing.
    #event_lbt_tmp = event_lbt_tmp[0:100]

    # I need to convert this to a list of numpy arrays.
    event_lbt = np.array([event_lbt_tmp[col].data for col in event_lbt_tmp.colnames])
    print(type(event_lbt))

    blend_rad = 1.0  # arcseconds
    blend_lbt = synthetic._calc_blends(bigpatch, c, event_lbt, blend_rad)

    # Fetch the first lens:
    lens_obj_id_1 = blend_lbt[0, 0]
    
    return
