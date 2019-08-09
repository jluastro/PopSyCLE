from popsycle import synthetic
from astropy.coordinates import SkyCoord # High-level coordinates
import astropy.coordinates as coord
from astropy.table import Table
import ebf
from astropy import units
import numpy as np
import h5py
import pdb
from multiprocessing import Pool
import itertools
import matplotlib.pyplot as plt
import time

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

def time_calc_event_time_loop():
    t1 = test_calc_event_time_loop(1)
    t4 = test_calc_event_time_loop(4)
    t8 = test_calc_event_time_loop(8)
    t12 = test_calc_event_time_loop(12)
    
    t_proc = np.array([t1, t4, t8, t12])
    t_thry = t_proc[0] * np.array([1, 1/4.0, 1/8.0, 1/12.0])
    n_proc = np.array([1, 4, 8, 12])
    
    fig = plt.figure(1, figsize=(6,6))
    plt.clf()
    plt.plot(n_proc, t_thry, '.', label = 'Ideal')
    plt.plot(n_proc, t_proc, '.', label = 'Reality')
    plt.xlabel('N processors')
    plt.ylabel('Time (sec)')
    plt.legend()
    plt.show()

def test_calc_event_time_loop(n_proc):
    t0 = time.time()
    
    test_data_dir = '/u/casey/scratch/work/microlens/galaxia_test/OGLE672/'
    hdf5_file = test_data_dir + 'OGLE672.h5'
    obs_time = 1000
    n_obs = 101
    radius_cut = 2
    theta_frac = 2
    blend_rad = 0.8

    # Initialize events_tmp and blends_tmp.
    events_tmp = None
    blends_tmp = None

    # Get the l and b from the HDF5 file.
    hf = h5py.File(hdf5_file, 'r')
    l_array = np.array(hf['long_bin_edges'])
    b_array = np.array(hf['lat_bin_edges'])
    hf.close()
    
    # Converts radius_cut from arcseconds into milliarcseconds
    radius_cut *= 1000.0

    # Set up the multiprocessing
    pool = Pool(n_proc)

    # Set up inputs to be able to be read by pool.map
    # Only do a few of the bins.
    nll = 7
    nbb = 7
    
    llbb = itertools.product(range(nll), range(nbb))

    reps = nll * nbb

    hd = itertools.repeat(hdf5_file, reps)
    ot = itertools.repeat(obs_time, reps)
    no = itertools.repeat(n_obs, reps)
    rc = itertools.repeat(radius_cut, reps)
    tf = itertools.repeat(theta_frac, reps) 
    br = itertools.repeat(blend_rad, reps)

    inputs = zip(llbb, hd, ot, no, rc, tf, br)

    ##########
    # Loop through galactic latitude and longitude bins. For each bin vertex, take
    # the nearest 4 bin samples and calculate microlensing events. We do this
    # to properly handle bin edges (i.e. a sliding window analysis of 2x2 bins).
    # Duplicate events are removed. 
    ##########
    # map_async? Make it wait till there are all done.
    results = pool.starmap(synthetic._calc_event_time_loop, inputs)

    pool.close()
    pool.join()

    results = [i for i in results if i is not None]
 
    # Is there a way for this to NOT involve a loop?????
    if len(results) > 0:
        for ii in range(len(results)):
            if events_tmp is not None:
                events_tmp = np.concatenate((events_tmp, results[ii][0]), axis=1)
            else:
                events_tmp = results[ii][0]

        for ii in range(len(results)):
            if blends_tmp is not None:
                blends_tmp = np.concatenate((blends_tmp, results[ii][1]), axis=1)
            else:
                blends_tmp = results[ii][1]

    t1 = time.time()

    runtime = t1 - t0

    return runtime
