import pytest

from popsycle import synthetic
from popsycle import ebf
from popsycle import binary_utils
import popsycle
from spisea.imf import multiplicity
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.coordinates as coord
from astropy.table import Table
from astropy import units
from astropy.io.misc import hdf5
import numpy as np
import h5py
from multiprocessing import Pool
import itertools
import matplotlib.pyplot as plt
import time
import filecmp
import os
from math import isclose
import resource
from copy import deepcopy
from scipy.spatial import cKDTree
from scipy.signal import find_peaks

masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcday = 1.0 * (3.086 * 10 ** 16) ** -1 * 86400.0
au_to_kpc = 4.848 * 10 ** -9

# Define some fixtures for each stage of synthetic analysis.
def galaxia():
    seed = 10

    output_root = 'data_test/test'

    test_filepath = os.path.dirname(__file__)
    galaxia_params = test_filepath + '/galaxyModelParams_PopSyCLEv3.txt'

    synthetic.write_galaxia_params(output_root=output_root,
                                   longitude=1.25,
                                   latitude=-2.65,
                                   area=0.001,
                                   seed=seed)

    synthetic.run_galaxia(output_root=output_root,
                          longitude=1.25,
                          latitude=-2.65,
                          area=0.001,
                          galaxia_galaxy_model_filename=galaxia_params,
                          seed=seed)

    return output_root

@pytest.fixture(name = 'galaxia', scope="module")
def galaxia_fixture():
    return galaxia()
    
def srun_galaxia(galaixa):
    # Srun and Mrun share the same galaxia file
    input_root = 'data_test/test'

    output_root = 'data_test/test_Srun'
    
    test_filepath = os.path.dirname(__file__)
    full_input_root = test_filepath + '/' + input_root
    full_output_root = test_filepath + '/' + output_root
    

    try:
        os.symlink(full_input_root + '_galaxia.log', full_output_root + '_galaxia.log')  
    except FileExistsError:
        pass
    try: 
        os.symlink(full_input_root + '_galaxia_params.txt', full_output_root + '_galaxia_params.txt')
    except FileExistsError:
        pass
    try:
        os.symlink(full_input_root + '.ebf', full_output_root + '.ebf')
    except FileExistsError:
        pass

    return output_root

@pytest.fixture(name = 'srun_galaxia', scope="module")
def srun_galaxia_fixture(galaxia):
    return srun_galaxia(galaxia)

def srun_popsyn(srun_galaxia):
    seed = 10

    ebf_file = srun_galaxia + '.ebf'
    output_root = srun_galaxia


    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='SukhboldN20',
                              overwrite=True,
                              seed=seed, n_proc=1)

    return output_root

@pytest.fixture(name = 'srun_popsyn', scope="module")
def srun_popsyn_fixture(srun_galaxia):
    return srun_popsyn(srun_galaxia)

def srun_calc_events(srun_popsyn):
    seed = 10

    hdf5_file = srun_popsyn + '.h5'
    output_root = srun_popsyn

    synthetic.calc_events(hdf5_file=hdf5_file,
                          output_root2=output_root,
                          radius_cut=2,
                          obs_time=1000,
                          n_obs=11,
                          theta_frac=2,
                          blend_rad=0.65,
                          overwrite=True,
                          n_proc=1)

    return output_root

@pytest.fixture(name = 'srun_calc_events', scope="module")
def srun_calc_events_fixture(srun_popsyn):
    return srun_calc_events(srun_popsyn)

def srun_refine_events(srun_calc_events):
    seed = 0

    input_root = srun_calc_events

    synthetic.refine_events(input_root=input_root,
                            filter_name='I',
                            photometric_system='ubv',
                            red_law='Damineli16',
                            overwrite=True,
                            output_file='default')

    output_root = input_root + '_refined_events_ubv_I_Damineli16'

    return output_root

@pytest.fixture(name = 'srun_refine_events', scope="module")
def srun_refine_events_fixture(srun_calc_events):
    return srun_refine_events(srun_calc_events)

def mrun_galaxia(galaixa):
    # Srun and Mrun share the same galaxia file
    input_root = 'data_test/test'
    output_root = 'data_test/test_Mrun'
    
    test_filepath = os.path.dirname(__file__)
    full_input_root = test_filepath + '/' + input_root
    full_output_root = test_filepath + '/' + output_root

    try:
        os.symlink(full_input_root + '_galaxia.log', full_output_root + '_galaxia.log')
    except FileExistsError:
        pass
    try:
        os.symlink(full_input_root + '_galaxia_params.txt', full_output_root + '_galaxia_params.txt')
    except FileExistsError:
        pass
    try:
        os.symlink(full_input_root + '.ebf', full_output_root + '.ebf')
    except FileExistsError:
        pass

    return output_root

@pytest.fixture(name = 'mrun_galaxia', scope="module")
def mrun_galaxia_fixture(galaxia):
    return mrun_galaxia(galaxia)

def mrun_popsyn(mrun_galaxia):
    seed = 10

    ebf_file = mrun_galaxia + '.ebf'
    output_root = mrun_galaxia

    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)

    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='SukhboldN20',
                              multiplicity=multi_obj,
                              overwrite=True,
                              seed=seed,
                              n_proc=1)

    return output_root

@pytest.fixture(name = 'mrun_popsyn', scope="module")
def mrun_popsyn_fixture(mrun_galaxia):
    return mrun_popsyn(mrun_galaxia)

def mrun_calc_events(mrun_popsyn):
    seed = 10

    hdf5_file = mrun_popsyn + '.h5'
    hdf5_comp_file = mrun_popsyn + '_companions.h5'
    output_root = mrun_popsyn

    synthetic.calc_events(hdf5_file=hdf5_file,
                          output_root2=output_root,
                          radius_cut=2,
                          obs_time=1000,
                          n_obs=11,
                          theta_frac=2,
                          blend_rad=0.65,
                          hdf5_file_comp=hdf5_comp_file,
                          overwrite=True,
                          n_proc=1)

    return output_root

@pytest.fixture(name = 'mrun_calc_events', scope="module")
def mrun_calc_events_fixture(mrun_popsyn):
    return mrun_calc_events(mrun_popsyn)

def mrun_refine_events(mrun_calc_events):
    seed = 10
    
    input_root = mrun_calc_events

    synthetic.refine_events(input_root=input_root,
                            filter_name='I',
                            photometric_system='ubv',
                            red_law='Damineli16',
                            hdf5_file_comp=input_root + '_companions.h5',
                            overwrite=True,
                            output_file='default',
                            seed = seed)

    output_root = input_root + '_refined_events_ubv_I_Damineli16'

    return output_root

@pytest.fixture(name = 'mrun_refine_events', scope="module")
def mrun_refine_events_fixture(mrun_calc_events):
    return mrun_refine_events(mrun_calc_events)

def mrun_refine_binary(mrun_refine_events):
    input_root = mrun_refine_events
    events_prim_fits = input_root + '.fits'
    events_comp_fits = input_root + '_companions.fits'
    synthetic.refine_binary_events(events_prim_fits,
                                   events_comp_fits,
                                   filter_name='I',
                                   photometric_system='ubv', overwrite=True,
                                   output_file='default', save_phot=False)

    output_root = input_root + '_rb'

    return output_root

@pytest.fixture(name = 'mrun_refine_binary', scope="module")
def mrun_refine_binary_fixture(mrun_refine_events):
    return mrun_refine_binary(mrun_refine_events)

def mrun_big_galaxia():
    seed = 10

    output_root = 'data_test/test_Mrun_big'

    test_filepath = os.path.dirname(__file__)
    galaxia_params = test_filepath + '/galaxyModelParams_PopSyCLEv3.txt'

    synthetic.write_galaxia_params(output_root=output_root,
                                   longitude=1.25,
                                   latitude=-2.65,
                                   area=0.0034,
                                   seed=seed)

    synthetic.run_galaxia(output_root=output_root,
                          longitude=1.25,
                          latitude=-2.65,
                          area=0.0034,
                          galaxia_galaxy_model_filename=galaxia_params,
                          seed=seed)

    return output_root

def mrun_big_popsyn():
    seed = 10

    output_root = 'data_test/test_Mrun_big'
    ebf_file = output_root + '.ebf'

    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)

    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='SukhboldN20',
                              multiplicity=multi_obj,
                              overwrite=True,
                              seed=seed, n_proc=4)

    max_mem_Mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1.0e6
    print(f'Max Memory Used: {max_mem_Mb} Mb')

    return output_root

def mrun_0p034deg_galaxia_NERSC():
    seed = 10

    output_root = 'data_test/test_Mrun_0p034deg'

    test_filepath = os.path.dirname(__file__)
    galaxia_params = test_filepath + '/galaxyModelParams_PopSyCLEv3.txt'

    synthetic.write_galaxia_params(output_root=output_root,
                                   longitude=1.25,
                                   latitude=-2.65,
                                   area=0.034,
                                   seed=seed)

    synthetic.run_galaxia(output_root=output_root,
                          longitude=1.25,
                          latitude=-2.65,
                          area=0.034,
                          galaxia_galaxy_model_filename=galaxia_params,
                          seed=seed)

    return output_root

def mrun_0p034deg_popsyn_NERSC():
    seed = 10

    output_root = 'data_test/test_Mrun_0p034deg'
    ebf_file = output_root + '.ebf'

    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)

    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/global/homes/n/nsabrams/uLens/code/src/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='SukhboldN20',
                              multiplicity=multi_obj,
                              overwrite=True,
                              seed=seed, n_proc=6)

    max_mem_Mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1.0e6
    print(f'Max Memory Used: {max_mem_Mb} Mb')

    return output_root

def remake_correct_data():
    import os, glob

    dir_data_test = os.path.dirname(__file__) + 'data_test/'
    dir_data_corr = os.path.dirname(__file__) + 'data_corr/'

    old_srun_files = glob.glob(dir_data_test + 'test_Srun*')
    old_mrun_files = glob.glob(dir_data_test + 'test_Mrun*')

    for ff in old_srun_files:
        os.copyfile(ff, dir_data_corr)

    for ff in old_mrun_files:
        os.copyfile(ff, dir_data_corr)

    return


def test_galactic_to_heliocentric():
    """
    Test the conversion from Galactic to Heliocentric positions.

    NOTE THAT HELIOCENTRIC COORDINATES ARE NOT THE SAME AS
    ASTROPY'S DEFAULT GALACTOCENTRIC COORDINATES!!!
    
    THEY ARE BOTH RECTANGULAR, BUT ONE IS CENTERED ON THE SUN,
    THE OTHER IS CENTERED ON THE GC!
    """
    
    # The arrays px, py, pz, glat, glon come from Galaxia.
    # They are just the first 10 entries of an EBF file it made.
    # For internal reference: the file is
    # /g2/scratch/casey/papers/microlens_2019/popsycle_rr_files/OGLE150_v2.ebf

    # Heliocentric 3-D positions (kpc)
    px = np.array([1.1423072 , 1.5105276 , 0.63417935, 0.7515631 , 1.1420373 ,
                   1.1334454 , 1.1531974 , 1.2154862 , 1.2258973 , 1.3694662 ])

    py = np.array([0.03450311, 0.04551123, 0.01913031, 0.02276252, 0.03423306,
                   0.0343769 , 0.03504916, 0.03649908, 0.03713397, 0.04156267])

    pz = np.array([-0.08887936, -0.11735878, -0.04948796, -0.05863307, -0.08901644,
                   -0.08859035, -0.08981982, -0.09497609, -0.09539149, -0.10696827])

    # Galactic 3-D positions (deg and kpc)
    glat = np.array([-4.4470243, -4.4406023, -4.4599943, -4.458851 , -4.454934 ,
                     -4.467121 , -4.451596 , -4.465916 , -4.4473953, -4.4642286])

    glon = np.array([1.7300799, 1.7257639, 1.7278297, 1.734782 , 1.716952 ,
                     1.7372237, 1.7408566, 1.7199831, 1.7350312, 1.7383676])

    rad = np.array([1.146279 , 1.5157632, 0.6363949, 0.7541903, 1.1460127,
                    1.1374218, 1.157221 , 1.2197374, 1.2301637, 1.374266 ])

    # Convert the Heliocentric positions to Galactic positions
    r, b, l = synthetic.heliocentric_to_galactic(px, py, pz)

    np.testing.assert_almost_equal(r, rad, 5)
    np.testing.assert_almost_equal(b, glat, 5)
    np.testing.assert_almost_equal(l, glon, 5)

    ##########
    # Astropy.
    # 
    # Note: This is not a test because we're showing
    # that Heliocentric is not the same as Galactocentric...
    # This is just for reference... so it's just comments.
    ##########
    # c_xyz = coord.SkyCoord(x=px * u.kpc,
    #                        y=py * u.kpc,
    #                        z=pz * u.kpc,
    #                        frame=coord.Galactocentric)
    # 
    # c_gal = c_xyz.transform_to(coord.Galactic)
    # 
    # print(np.min(c_gal.b.value - glat))
    # print(np.min(c_gal.l.value - glon))
    # print(np.min(c_gal.distance.value - rad))
    
def test_galactic_to_heliocentric_1():
    """
    Test that the formula gives the right things in 3 simple geometry cases.
    """
    r1, b1, l1 = synthetic.heliocentric_to_galactic(1, 0, 0)
    np.testing.assert_equal([r1, b1, l1], [1, 0, 0])
    
def test_galactic_to_heliocentric_2():
    """
    Test that the formula gives the right things in 3 simple geometry cases.
    """
    r2, b2, l2 = synthetic.heliocentric_to_galactic(0, 1, 0)
    np.testing.assert_equal([r2, b2, l2], [1, 0, 90])

def test_galactic_to_heliocentric_3():
    """
    Test that the formula gives the right things in 3 simple geometry cases.
    """
    r3, b3, l3 = synthetic.heliocentric_to_galactic(0, 0, 1)
    np.testing.assert_equal([r3, b3, l3], [1, 90, 0])



@pytest.mark.xfail
def test_calc_blends(srun_popsyn, srun_calc_events):
    hdf5_file = srun_popsyn + '.h5'
    events_file = srun_calc_events

    test_h5 = h5py.File(hdf5_file, 'r')
    subfield_list = list(test_h5.keys())[1:-2]

    # just work on the first field.
    bigpatch1 = hdf5.read_table_hdf5(test_h5[subfield_list[0]])
    bigpatch2 = hdf5.read_table_hdf5(test_h5[subfield_list[0]])
    bigpatch3 = hdf5.read_table_hdf5(test_h5[subfield_list[0]])
    bigpatch4 = hdf5.read_table_hdf5(test_h5[subfield_list[0]])
    bigpatch = np.concatenate([bigpatch1, bigpatch2, bigpatch3, bigpatch4])

    time = 0.0
    r_t = bigpatch['rad'] + time * bigpatch['vr'] * kms_to_kpcday  # kpc
    b_t = bigpatch['glat'] + time * bigpatch['mu_b'] * masyr_to_degday  # deg
    l_t = bigpatch['glon'] + time * (bigpatch['mu_lcosb'] / np.cos(np.radians(bigpatch['glat']))) * masyr_to_degday  # deg

    ##########
    # Define coordinates.
    ##########
    c = SkyCoord(frame='galactic', l=l_t * units.deg, b=b_t * units.deg)
    foo = coord.match_coordinates_sky(c, c, nthneighbor=2)
    print(foo)

    # Load up the events table
    event_lbt_tmp = Table.read(events_file)

    # Trim it down for faster processing.
    #event_lbt_tmp = event_lbt_tmp[0:100]

    # Convert this to a list of numpy arrays via pandas.
    event_lbt_tmp_df = event_lbt_tmp.to_pandas()
    event_lbt = event_lbt_tmp_df.to_records()

    blend_rad = 20.0  # arcseconds
    blend_lbt = synthetic._calc_blends(bigpatch, c, event_lbt, blend_rad)

    # Fetch the first lens:
    lens_obj_id_1 = blend_lbt[0, 0]

    return

# FIXME
@pytest.mark.xfail
def test_calc_blends_old():
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
    #bigpatch = np.hstack((hf[name00], hf[name01], hf[name10], hf[name11]))
    hf.close()
    print(bigpatch.shape)

    time = 0.0
    r_t = bigpatch[9] + time * bigpatch[12] * kms_to_kpcday  # kpc
    b_t = bigpatch[10] + time * bigpatch[13] * masyr_to_degday  # deg
    l_t = bigpatch[11] + time * (bigpatch[14] / np.cos(np.radians(bigpatch[10]))) * masyr_to_degday  # deg

    ##########
    # Define coordinates.
    ##########
    c = SkyCoord(frame='galactic', l=l_t * units.deg, b=b_t * units.deg)
    foo = coord.match_coordinates_sky(c, c, nthneighbor=2)

    # Load up the events table
    event_lbt_tmp = Table.read(event_file)

    # Trim it down for faster processing.
    # event_lbt_tmp = event_lbt_tmp[0:100]

    # Convert this to a list of numpy arrays via pandas.
    event_lbt_tmp_df = event_lbt_tmp.to_pandas()
    event_lbt = event_lbt_tmp_df.to_records()
    
    #event_lbt_tmp_types = [np.dtype(event_lbt_tmp[col].data[0]) for col in event_lbt_tmp.colnames]
    #print(list(tuple(zip(event_lbt_tmp.colnames, event_lbt_tmp_types))))
    #print([event_lbt_tmp[col].data for col in event_lbt_tmp.colnames])
    #event_lbt = np.array([event_lbt_tmp[col].data for col in event_lbt_tmp.colnames], dtype = list(tuple(zip(event_lbt_tmp.colnames, event_lbt_tmp_types))))
    #print(event_lbt_tmp['glat_L'])
    #print(event_lbt['glat_L'])

    #print(type(event_lbt))
    print(bigpatch[0])

    blend_rad = 1.0  # arcseconds
    blend_lbt = synthetic._calc_blends(bigpatch, c, event_lbt, blend_rad)

    # Fetch the first lens:
    lens_obj_id_1 = blend_lbt[0, 0]

    return

# FIXME
# def time_calc_event_time_loop():
#     t1 = test_calc_event_time_loop(1)
#     t4 = test_calc_event_time_loop(4)
#     t8 = test_calc_event_time_loop(8)
#     t12 = test_calc_event_time_loop(12)
#
#     t_proc = np.array([t1, t4, t8, t12])
#     t_thry = t_proc[0] * np.array([1, 1 / 4.0, 1 / 8.0, 1 / 12.0])
#     n_proc = np.array([1, 4, 8, 12])
#
#     fig = plt.figure(1, figsize=(6, 6))
#     plt.clf()
#     plt.plot(n_proc, t_thry, '.', label='Ideal')
#     plt.plot(n_proc, t_proc, '.', label='Reality')
#     plt.xlabel('N processors')
#     plt.ylabel('Time (sec)')
#     plt.legend()
#     plt.show()

# FIXME
def calc_event_time_loop(n_proc):
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

    #    # Set up inputs to be able to be read by pool.map
    #    # Only do a few of the bins.
    nll = range(10, 18, 1)
    nbb = range(10, 11, 1)

    llbb = itertools.product(nll, nbb)

    reps = len(nll) * len(nbb)

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
                events_tmp = np.concatenate((events_tmp, results[ii][0]),
                                            axis=1)
            else:
                events_tmp = results[ii][0]

        for ii in range(len(results)):
            if blends_tmp is not None:
                blends_tmp = np.concatenate((blends_tmp, results[ii][1]),
                                            axis=1)
            else:
                blends_tmp = results[ii][1]

    t1 = time.time()

    runtime = t1 - t0

    return runtime


@pytest.mark.xfail
def test_calc_events():
    test_data_dir = '/u/casey/scratch/work/microlens/galaxia_test/OGLE672/'
    hdf5_file = test_data_dir + 'OGLE672.h5'
    output_root2 = 'OGLE672_TEST'
    obs_time = 1000
    n_obs = 101
    radius_cut = 2
    theta_frac = 2
    blend_rad = 0.8
    n_proc = 1
    overwrite = False

    t0 = time.time()

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
    # Should I use starmap_async? 
    results = pool.starmap(synthetic._calc_event_time_loop, inputs)

    pool.close()
    pool.join()

    # Remove all the None values
    # (occurs for patches with less than 10 objects)
    results = [i for i in results if i is not None]

    # Remove all the None values
    # (occurs for patches with less than 10 objects)
    results = [i for i in results if i is not None]

    results_ev = []
    results_bl = []

    for ii in range(len(results)):
        if results[ii] is not None:
            if results[ii][0] is not None:
                results_ev.append(results[ii][0])
            if results[ii][1] is not None:
                results_bl.append(results[ii][1])

    events_tmp = np.concatenate(results_ev, axis=1)
    blends_tmp = np.concatenate(results_bl, axis=1)

    # Convert the events numpy array into an Astropy Table for easier consumption.
    # The dimensions of events_final_table are 52 x Nevents
    if events_tmp is not None:
        events_tmp = synthetic.unique_events(events_tmp)
        events_final = Table(events_tmp.T,
                             names=('zams_mass_L', 'rem_id_L', 'mass_L',
                                    'px_L', 'py_L', 'pz_L',
                                    'vx_L', 'vy_L', 'vz_L',
                                    'rad_L', 'glat_L', 'glon_L',
                                    'vr_L', 'mu_b_L', 'mu_lcosb_L',
                                    'age_L', 'popid_L', 'ubv_k_L', 'ubv_i_L',
                                    'exbv_L', 'obj_id_L',
                                    'ubv_j_L', 'ubv_u_L', 'ubv_r_L',
                                    'ubv_b_L', 'ubv_h_L', 'ubv_v_L',
                                    'zams_mass_S', 'rem_id_S', 'mass_S',
                                    'px_S', 'py_S', 'pz_S',
                                    'vx_S', 'vy_S', 'vz_S',
                                    'rad_S', 'glat_S', 'glon_S',
                                    'vr_S', 'mu_b_S', 'mu_lcosb_S',
                                    'age_S', 'popid_S', 'ubv_k_S', 'ubv_i_S',
                                    'exbv_S', 'obj_id_S',
                                    'ubv_j_S', 'ubv_u_S', 'ubv_r_S',
                                    'ubv_b_S', 'ubv_h_S', 'ubv_v_S',
                                    'theta_E', 'u0', 'mu_rel', 't0'))

    if blends_tmp is not None:
        blends_tmp = synthetic.unique_blends(blends_tmp)
        blends_final = Table(blends_tmp.T, names=('obj_id_L', 'obj_id_S',
                                                  'zams_mass_N', 'rem_id_N',
                                                  'mass_N',
                                                  'px_N', 'py_N', 'pz_N',
                                                  'vx_N', 'vy_N', 'vz_N',
                                                  'rad_N', 'glat_N', 'glon_N',
                                                  'vr_N', 'mu_b_N',
                                                  'mu_lcosb_N',
                                                  'age_N', 'popid_N',
                                                  'ubv_k_N', 'ubv_i_N',
                                                  'exbv_N', 'obj_id_N',
                                                  'ubv_j_N', 'ubv_u_N',
                                                  'ubv_r_N',
                                                  'ubv_b_N', 'ubv_h_N',
                                                  'ubv_v_N',
                                                  'sep_LN'))

    if events_tmp is None:
        print('No events!')
        return

    t1 = time.time()

    # Save out file 
    events_final.write(output_root2 + '_events.fits', overwrite=overwrite)
    blends_final.write(output_root2 + '_blends.fits', overwrite=overwrite)

    print('Total runtime: {0:f} s'.format(t1 - t0))

    return

def test_add_multiples():
    """
    Checks that _add_multiples() is making the companions point to the correct Galaxia objects
    and eliminates the correct elements.
    """
    # Fake version of Galaxia star masses
    star_masses_check = np.linspace(0, 5, 15).tolist()
    
    class FakeCluster(object):
        """
        Makes a "fake" cluster which is the necessary parts
        of the star_systems and companions tables
        """
        
        # Makes mass slightly off from the star_zams_masses
        mass = np.linspace(0.1, 5.1, 10)
        
        # All odd indices are Multiple systems
        isMulti = np.empty((len(mass),))
        isMulti[::2] = False
        isMulti[1::2] = True
        
        # Makes first two white dwarfs and the rest stars
        phase = np.empty((len(mass),))
        phase[0:2] = 101
        phase[2:] = 1.0

        # Set mass and zams_mass equal for this fake cluster.
        star_systems = Table([mass, mass, isMulti, phase],
                             names=['mass', 'zams_mass', 'isMultiple', 'phase'])

        # Make companions.
        system_idx = [1,1,1,3,5,5,7,9]
        comp_phase = np.ones((len(system_idx)), dtype=int)
        comp_mass = np.random.rand(len(system_idx)) * mass[system_idx]  # uniform mass ratio.
        companions = Table([system_idx, comp_phase, comp_mass, comp_mass],
                             names=['system_idx', 'phase', 'mass', 'zams_mass'])
        
    cluster_check = FakeCluster()
        
    companion_check = synthetic._add_multiples(star_masses_check, cluster_check)
        
    # 1,1,1 should be eliminated because the primary is a WD
    # 9 should be eliminated because it's too big
    # 3, 5, 5, and 7 should be pointed to their star_zams_masses counterpart
    system_idx_correct = [5,9,9,12]

    if not np.array_equal(system_idx_correct, companion_check['system_idx']):
        raise Exception("_add_multiples() is behaving unexpectedly")
    
    return

def test_system_mass(mrun_popsyn):
    """
    Test whether the system mass (current) is really the sum of the individual component masses.

    Run working_test_all_Mrun() first.
    """

    test_filepath = os.path.dirname(__file__)
    seed = 10

    output_root = mrun_popsyn
    
    test_h5 = h5py.File(output_root + '.h5', 'r')
    test_comp_h5 = h5py.File(output_root + '_companions.h5', 'r')

    subfield_list = list(test_h5.keys())[1:-2]
    
    for ff in range(len(subfield_list)):
    #for ff in range(0, 1):
        field = subfield_list[ff]

        # Convert to astropy tables for easier work. These are small enough.
        tab_syst = hdf5.read_table_hdf5(test_h5[field])
        tab_comp = hdf5.read_table_hdf5(test_comp_h5[field])

        # Make sure there is a 'mass' and 'zams_mass' column
        # in both primary and companion tables.
        assert 'mass' in tab_syst.colnames
        assert 'mass' in tab_comp.colnames
        assert 'zams_mass' in tab_syst.colnames
        assert 'zams_mass' in tab_comp.colnames

        # Loop through each object in the systems table and find
        # the associated companions. Do it the slow way to be absolutely sure.
        n_comp_checked = 0
        #for ss in range(len(tab_syst)):
        if len(tab_syst) == 0:
            continue
            
        for ss in range(1000):
            primary = tab_syst[ss]

            if primary['isMultiple']:
                cdx = np.where(tab_comp['system_idx'] == primary['obj_id'])[0]
                companions = tab_comp[cdx]

                # Sum up the primary + companion current masses.
                mass = np.sum(companions['mass'])
                mass += primary['mass']

                # Sum up the primary + companion zams masses.
                zmass = np.sum(companions['zams_mass'])
                zmass += primary['zams_mass']

                # While we are here, check the current mass is < the zams mass.
                # Adding + 10**-2 since at very low masses sometimes isochrone
                # puts zams_mass < mass
                mass_epsilon = 1e-2
                assert primary['mass'] <= primary['zams_mass'] + mass_epsilon
                for comp in companions:
                    assert comp['mass'] <= comp['zams_mass'] + mass_epsilon
                    
                # check for companion < primary mass
                mass_ratio_epsilon = 1.1   # comp mass is allowed to be 10% greater than primary
                for comp in companions:
                    assert comp['zams_mass'] <= primary['zams_mass'] * mass_ratio_epsilon

                # Confirm that the sum of masses equals the system mass.
                assert isclose(mass, primary['systemMass'], abs_tol=1e-5)

                # Confirm that the current system mass is less than the zams system mass.
                assert mass <= zmass + mass_epsilon

                # Keep track of how many we checked.
                n_comp_checked += 1

        print(f'Field {field}: Successfully checked system masses of {n_comp_checked} multiple systems.')
        
    test_h5.close()
    test_comp_h5.close()
    
    return


def test_system_luminosity(mrun_galaxia):
    """
    Checks that the system luminosity is the
    primary + companions with compact objects not having luminosity.
    Not integration test since primary luminosities are not kept.
    """
    
    primary_star_dict, co_dict, star_dict, companions_table = stripped_perform_pop_syn_for_system_luminosity_test(mrun_galaxia)
    
    #Group companions by system_idx
    grouped_companions = companions_table.group_by(['system_idx'])
    companions_system_m_ubv_I = grouped_companions['m_ubv_I'].groups.aggregate(binary_utils.add_magnitudes)
    grouped_system_idxs = np.array(grouped_companions.groups.keys['system_idx'])
    
    # Returns the intersecting obj_ids/system_idxs, indices in co_table or star_dict that correspond with the overlap, 
    # and indices in companion_table that overlap (this last one should be just np.arange(len(companion_table)))
    # DOES NOT TEST WDs since some are luminous and I have no reference magnitudes for them!!!
    NS_BH_idxs = np.where(co_dict['rem_id'] >= 102)
    co_idx_w_companions = np.intersect1d(np.array(co_dict['obj_id'][NS_BH_idxs]), grouped_system_idxs, return_indices = True)
    star_idx_w_companions = np.intersect1d(np.array(star_dict['obj_id']), grouped_system_idxs, return_indices = True)
    
    # checks the matching was done properly
    assert np.array_equiv(star_dict['obj_id'][star_idx_w_companions[1]], grouped_system_idxs[star_idx_w_companions[2]])
    assert np.array_equiv(co_dict['obj_id'][NS_BH_idxs][co_idx_w_companions[1]], grouped_system_idxs[co_idx_w_companions[2]])
    
    calc_primary_luminosity = binary_utils.primary_mag_from_system_mag(star_dict['ubv_I'][star_idx_w_companions[1]],
                                                                       companions_system_m_ubv_I[star_idx_w_companions[2]])
    
    # checks the primary luminosity is the same as 
    # the system luminosity - companion luminosity
    assert np.allclose(calc_primary_luminosity, primary_star_dict['ubv_I'][star_idx_w_companions[1]])
    
    calc_primary_luminosity_co = binary_utils.primary_mag_from_system_mag(co_dict['ubv_I'][NS_BH_idxs][co_idx_w_companions[1]],
                                                                       companions_system_m_ubv_I[co_idx_w_companions[2]])
    

    # checks that all the compact object 
    # primary luminosities are nan
    assert np.all(np.isnan(calc_primary_luminosity_co))
    
    return

def stripped_perform_pop_syn_for_system_luminosity_test(mrun_galaxia):
    """
    Mock version of perform_pop_syn that will also return the 
    primary magnitudes. Runs for an arbitrary pop_id, age bin, and
    metallicity bin.
    If the way in which the bins are chosen are changed, this function
    must also be changed.
    """
    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)
    
    test_filepath = os.path.dirname(__file__)
    seed = 1
    
    pid = 7
    age_min = 10.041392
    age_max = 10.041393
    feh_min = -1.500
    feh_max = 0.00
    
    additional_photometric_systems = None

    ebf_file = mrun_galaxia + '.ebf'
    
    popid_array = ebf.read(ebf_file, '/popid')
    popid_idx = np.where(popid_array == pid)[0]
    age_array = ebf.read(ebf_file, '/age')
    age_idx = np.where((age_array[popid_idx] >= age_min) &
                               (age_array[popid_idx] < age_max))[0]
    metallicity_array = ebf.read(ebf_file, '/feh')
    feh_idx = np.where((metallicity_array[popid_idx[age_idx]] >= feh_min) &
                                   (metallicity_array[popid_idx[age_idx]] < feh_max))[0]
    bin_idx = popid_idx[age_idx[feh_idx]]
    
    num_stars_in_bin = len(bin_idx)
    
    exbv_arr4kdt, kdt_star_p = synthetic._make_extinction_kdtree(ebf_file, popid_idx[age_idx])
        
    star_dict = {}
    star_dict['age'] = age_array[bin_idx]
    star_dict['popid'] = popid_array[bin_idx]
    star_dict['obj_id'] = np.arange(len(bin_idx))
    
    synthetic._load_galaxia_into_star_dict(star_dict,
                                 bin_idx,
                                 ebf_file,
                                 additional_photometric_systems)
    
    primary_star_dict = {}
    primary_star_dict['ubv_J'] = deepcopy(star_dict['ubv_J'])
    primary_star_dict['ubv_H'] = deepcopy(star_dict['ubv_H'])
    primary_star_dict['ubv_K'] = deepcopy(star_dict['ubv_K'])
    primary_star_dict['ubv_U'] = deepcopy(star_dict['ubv_U'])
    primary_star_dict['ubv_I'] = deepcopy(star_dict['ubv_I'])
    primary_star_dict['ubv_B'] = deepcopy(star_dict['ubv_B'])
    primary_star_dict['ubv_V'] = deepcopy(star_dict['ubv_V'])
    primary_star_dict['ubv_R'] = deepcopy(star_dict['ubv_R'])
    
    mass_in_bin = np.sum(star_dict['mass'])
    
    stars_in_bin = {}
    for key, val in star_dict.items():
        stars_in_bin[key] = val

    from multiprocessing import Lock, Value
    next_id_stars_val_in = Value('i', 0)
    next_id_co_val_in = Value('i', num_stars_in_bin)
    mp_lock = Lock()
    synthetic._mp_init_worker(mp_lock, next_id_stars_val_in, next_id_co_val_in)

    cluster, _, _ = synthetic._make_cluster(iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                            log_age = np.mean([age_min, age_max]), 
                            currentClusterMass = mass_in_bin, 
                            multiplicity=multi_obj, 
                            IFMR = 'SukhboldN20', 
                            feh = np.mean([feh_min, feh_max]), 
                            seed=seed)
    
    co_dict, _ = synthetic._make_co_dict(np.mean([age_min, age_max]),
                                     cluster,
                                     stars_in_bin,
                                     kdt_star_p, exbv_arr4kdt,
                                     BH_kick_speed_mean=100,
                                     NS_kick_speed_mean=350,
                                     multiplicity=multi_obj,
                                     seed=seed)
    
    star_dict, companions_table = synthetic._make_companions_table(cluster=cluster,
                                                         star_dict=star_dict,
                                                         co_dict=co_dict)
    
    
    return primary_star_dict, co_dict, star_dict, companions_table

def test_add_mags():
    """
    Tests binary_utils.add_magnitudes() for a variety of cases both
    individual and for a mock table of companions and primaries.
    """
    # array of 3 magnitudes adds up to the correct value
    assert binary_utils.add_magnitudes([10, 10, 10]) == 8.807196863200843 
    
    # including one nan gives the sum not including the nan
    assert binary_utils.add_magnitudes([10, np.nan]) == 10
    
    # sum of all nans gives nan
    assert np.isnan(binary_utils.add_magnitudes([np.nan, np.nan])) == True
    
    # check a mock companion table grouped by system idx
    companions_table = Table([[0, 0, 1, 2, 2, 3, 3], 
                        [np.nan, 10, 11, np.nan, np.nan, 12, 13]], 
                       names=('system_idx', 'mags'))
    
    grouped_companions = companions_table.group_by(['system_idx'])
    companions_system_m = grouped_companions['mags'].groups.aggregate(binary_utils.add_magnitudes)
    
    assert companions_system_m[0] == 10
    assert companions_system_m[1] == 11
    assert np.isnan(companions_system_m[2]) == True
    assert companions_system_m[3] == 11.636148842226767
    
    return

def test_subtract_mags():
    """
    Tests binary_utils.subtract_magnitudes() for a variety of cases both
    individual and for a mock table of companions and primaries.
    If primary is dimmer than the companion, it will give a warning and a
    nan magnitude.
    """
    
    # brighter star - dimmer star gives the correct value
    assert binary_utils.subtract_magnitudes(9, 10) == 9.551202076354773
    
    # magnitude - nan = magnitude
    assert binary_utils.subtract_magnitudes(9, np.nan) == 9.0
    
    # magnitude - magnitude = nan
    assert np.isnan(binary_utils.subtract_magnitudes(9, 9)) == True
    
    # nan - magnitude = nan
    assert np.isnan(binary_utils.subtract_magnitudes(np.nan, 9)) == True
    
    # nan - nan = nan
    assert np.isnan(binary_utils.subtract_magnitudes(np.nan, np.nan)) == True
    
    companions_table = Table([[0, 0, 1, 2, 2, 3, 3], 
                        [np.nan, 10, 11, np.nan, np.nan, 12, 13]], 
                       names=('system_idx', 'mags'))
    test_table_prims = np.array([np.nan, 9, np.nan, 10])
    
    grouped_companions = companions_table.group_by(['system_idx'])
    companions_system_m = grouped_companions['mags'].groups.aggregate(binary_utils.add_magnitudes)
    grouped_system_idxs = np.array(grouped_companions.groups.keys['system_idx'])
    
    prims = test_table_prims[grouped_system_idxs]
    
    subtracted_mags = binary_utils.subtract_magnitudes(prims, companions_system_m)
    
    assert np.isnan(subtracted_mags[0]) == True
    assert subtracted_mags[1] == 9.187350918581537
    assert np.isnan(subtracted_mags[2]) == True
    assert subtracted_mags[3] == 10.271972084483362
    
    return

def test_no_nan_companions(mrun_popsyn):
    """
    Tests that there are no nan companions that weren't
    assigned isochrones (happens for low mass companions often
    with high metalicity + old)

    Checks systemMass and mass are not nan and that N_companions
    is correct (since it was modified)
    """
    
    output_root = mrun_popsyn
    test_h5 = h5py.File(output_root + '.h5', 'r')
    test_comp_h5 = h5py.File(output_root + '_companions.h5', 'r')
    
    subfield_list = list(test_h5.keys())[1:-2]
    
    for ff, subfield in enumerate(subfield_list):
        systemMasses = test_h5[subfield]['systemMass']
        companion_masses = test_comp_h5[subfield]['mass']
        
        assert sum(np.isnan(systemMasses)) == 0
        assert sum(np.isnan(companion_masses)) == 0
    
        del systemMasses, companion_masses
    
        # Checks number of companions matches for 1000 cases
        # since this was edited
        N_companions_all = test_h5[subfield]['N_companions']
        obj_id_all = test_h5[subfield]['obj_id']
        system_id_all = test_comp_h5[subfield]['system_idx']
        for ss in range(1000):
            N_companions = N_companions_all[ss]
            obj_id = obj_id_all[ss]
            companions = np.where(system_id_all == obj_id)[0]
    
            assert(N_companions == len(companions))
    
        del N_companions_all, obj_id_all, system_id_all
        


def test_single_CO_frac(srun_popsyn):
    """
    Checks that the CO fraction of objects greater than 0.1 Msun
    is about 8.2%
    """
    test_hdf5 = h5py.File(srun_popsyn + '.h5', 'r')
    lower_mass_cutoff = 0.1 #Msun
    CO, total, CO_frac = calc_CO_frac_mass_cutoff(test_hdf5, lower_mass_cutoff)
    
    precalc_CO_frac = 0.0822
    precalc_CO_number = 2071
    precalc_total_number = 25184
    precalc_error = precalc_CO_frac*np.sqrt((np.sqrt(precalc_CO_number)/precalc_CO_number)**2 + (np.sqrt(precalc_total_number)/precalc_total_number)**2)
    
    assert(np.abs(CO_frac - precalc_CO_frac) < precalc_error)
    
    test_hdf5.close()
    
    return


def calc_CO_frac_mass_cutoff(hdf5_file, lower_mass_cutoff):
    subfield_list = list(hdf5_file.keys())[1:-2]
    CO = 0
    total = 0
    for field in subfield_list:
        array = hdf5_file[field]
        CO += len(np.where((array['rem_id'] > 100) & (array['mass'] > lower_mass_cutoff))[0])
        total += len(np.where(array['mass'] > lower_mass_cutoff)[0])
        del array
    CO_frac = CO/total
    return  CO, total, CO_frac


def test_multiplicity_properties(mrun_popsyn):
    """
    Checks that the multiplicity fraction of objects > 0.5 Msun is about 47%
    and that the minimum semimajor axis is greater than 10^-2
    """
    test_hdf5 = h5py.File(mrun_popsyn + '.h5', 'r')
    lower_mass_cutoff = 0.5 #Msun
    multiplicity_frac, multiples, total = calc_multiplicity_frac_mass_cutoff(test_hdf5, lower_mass_cutoff)
    
    precalc_mult_frac = 0.4727
    precalc_mult_number = 2052
    precalc_total_number = 4341
    precalc_error = precalc_mult_frac*np.sqrt((np.sqrt(precalc_mult_number)/precalc_mult_number)**2 + (np.sqrt(precalc_total_number)/precalc_total_number)**2)
    
    assert(np.abs(multiplicity_frac - precalc_mult_frac) < precalc_error)
    
    test_hdf5.close()

    test_hdf5_comp = h5py.File(mrun_popsyn + '_companions.h5', 'r')
    min_log_semimajor_axis = calc_min_semimajor_axis(test_hdf5_comp)
    
    min_log_semimajor_axis_imposed = -2
    
    assert(min_log_semimajor_axis_imposed < min_log_semimajor_axis)

    test_hdf5_comp.close()
    
    return

def calc_multiplicity_frac_mass_cutoff(hdf5_file, lower_mass_cutoff):
    subfield_list = list(hdf5_file.keys())[1:-2]
    multiples = 0
    total = 0
    for field in subfield_list:
        array = hdf5_file[field]
        multiples += len(np.where((array['isMultiple'] == 1) &
                                  (array['mass'] > lower_mass_cutoff))[0])
        total += len(np.where(array['mass'] > lower_mass_cutoff)[0])
        del array
    multiple_frac = multiples/total
    return  multiple_frac, multiples, total

def calc_min_semimajor_axis(hdf5_file):
    subfield_list = list(hdf5_file.keys())[1:-2]
    min_semimajor_axis = np.nan
    for field in subfield_list:
        array = hdf5_file[field]
        array_min_semimajor_axis = min(array['log_a'])
        if np.isnan(min_semimajor_axis):
            min_semimajor_axis = array_min_semimajor_axis
        elif array_min_semimajor_axis < min_semimajor_axis:
            min_semimajor_axis = array_min_semimajor_axis
            
        del array
    
    return  min_semimajor_axis




def test_binary_angles(mrun_refine_events):
    test_companions_table = Table.read(mrun_refine_events + '_companions.fits')
    
    alphas = test_companions_table['alpha']
    phis = test_companions_table['phi']
    phi_pi_Es = test_companions_table['phi_pi_E']
    
    assert(min(alphas) >= 0)
    assert(max(alphas) <= 360)
    assert(min(phis) >= 0)
    assert(max(phis) <= 360)
    assert(min(phi_pi_Es) >= 0)
    assert(max(phi_pi_Es) <= 360)
    
    return

# FIXME
@pytest.mark.xfail
def test_refine_binary_events_multiple_lightcurves(mrun_refine_binary):
    """
    Makes sure refine binary events properly 
    makes the right number of lightcurves for event (i.e. 2 lightcurves
    for PSBL and 4 lightcurves for BSBL). Also that it picks the correct
    one as the primary one (highest delta m).
    """
    
    test_events = Table.read(mrun_refine_binary + '.fits')

def test_refine_binary_events_psbl_lightcurve():
    """
    Generates a lightcurve that should have 4 peaks.
    Check by eye against 'tests/data_correct/psbl_4_peaks_example.png'
    """
    
    psbl_param_dict = {'raL': 265,
                       'decL' : -30,
                       'mL1' : 4,
                       'mL2' : 1,
                       't0' : 0,
                       'xS0': np.array([0, 0]),
                       'beta' : 0.001,
                       'muL' : np.array([0,0]),
                       'muS' : np.array([0,1]),
                       'dL' : 4000,
                       'dS' : 8000,
                       'sep' : 1,
                       'alpha' : 180,
                       'mag_src' : 22,
                      'b_sff' : 1}
    
    psbl = synthetic.psbl_model_gen(psbl_param_dict)
    
    duration=1000 # days
    time_steps=5000
    tmin = psbl.t0 - (duration / 2.0)
    tmax = psbl.t0 + (duration / 2.0)
    dt = np.linspace(tmin, tmax, time_steps)
    img, amp = psbl.get_all_arrays(dt)
    phot = psbl.get_photometry(dt, amp_arr=amp)
    
    peaks, _ = find_peaks(-phot, prominence = 10e-5, width =1) 
    
    assert(len(peaks) == 4)
    
    plt.plot(dt, phot.data)
    plt.plot(dt[peaks], phot.data[peaks], marker = '.', linestyle = 'None')
    plt.gca().invert_yaxis()
    plt.savefig(popsycle.__path__[0] + '/tests/data_test/psbl_4_peaks_example.png')
    
    return

def test_bspl_single_luminous_source_one_peak():
    """
    Makes sure that a BSPL event with one luminous source
    and one dark source only has a single peak (picked parameters
    such that not high parallax)
    """
    
    bspl_param_dict = {'raL': 265,
                       'decL' : -30,
                       'mL' : 4,
                       't0' : 0,
                       'xS0': np.array([0, 0]),
                       'beta' : 0.001,
                       'muL_E' : 0,
                       'muL_N' : 0,
                       'muS_E' : 0,
                       'muS_N' : 1,
                       'dL' : 4000,
                       'dL_dS' : 0.5,
                       'sep' : 1,
                       'alpha' : 270,
                       'mag_src_pri' : 22,
                         'mag_src_sec' : np.nan,
                      'b_sff' : 1}
    bspl = synthetic.bspl_model_gen(bspl_param_dict)
    
    duration=1000 # days
    time_steps=5000
    tmin = bspl.t0 - (duration / 2.0)
    tmax = bspl.t0 + (duration / 2.0)
    dt = np.linspace(tmin, tmax, time_steps)
    phot = bspl.get_photometry(dt)
    
    peaks, _ = find_peaks(-phot, prominence = 10e-5, width =1) 
    
    assert(len(peaks) == 1)
    
    return

#FIXME
@pytest.mark.xfail
def test_psbl_parallax():
    """
    Test psbl lightcurve model to make sure in a 
    bulge field we only have E component of parallax.
    This is to verify that l and b haven't been mixed up.
    """
    
    switch_l_and_b_psbl_parameter_dict = {'raL': 255.87254973655152,
                                          'decL': -26.747150542896915,
                                          'mL1': 0.5779610276222229,
                                          'mL2': 0.5368390723932487,
                                          't0': -251.56921080815064,
                                          'xS0': np.array([0, 0]),
                                          'beta': 0.07288389063108912,
                                          'muL': np.array([-2.93909249, -0.21650692]),
                                          'muS': np.array([-5.13010866, -6.9610909 ]),
                                          'dL': 1798.2280471728398,
                                          'dS': 7807.97229570119,
                                          'sep': 3.0858608224470054,
                                          'alpha': 348.88370965770963,
                                          'mag_src': 24.11974436088592,
                                          'b_sff': 0.006463304120418371}
    
    psbl_parameter_dict = {'raL': 274.61067095602954,
                              'decL': -22.847924243752118,
                              'mL1': 0.5779610276222229,
                              'mL2': 0.5368390723932487,
                              't0': -251.56921080815064,
                              'xS0': np.array([0, 0]),
                              'beta': 0.07288389063108912,
                              'muL': np.array([-2.87453565, -0.64975743]),
                              'muS': np.array([-4.04164504, -7.64459984]),
                              'dL': 1798.2280471728398,
                              'dS': 7807.97229570119,
                              'sep': 3.0858608224470054,
                              'alpha': 348.88370965770963,
                              'mag_src': 24.11974436088592,
                              'b_sff': 0.006463304120418371}
    
    model = synthetic.psbl_model_gen(test_dict)#psbl_parameter_dict)
    
    piE_E, piE_N = model.piE
    
    import pdb
    pdb.set_trace()
    
    return

