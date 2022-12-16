import pytest

from popsycle import synthetic
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

masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcday = 1.0 * (3.086 * 10 ** 16) ** -1 * 86400.0
au_to_kpc = 4.848 * 10 ** -9

# Define some fixtures for each stage of synthetic analysis.
@pytest.fixture
def srun_galaxia():
    seed = 1

    output_root = 'data_test/test_Srun'

    test_filepath = os.path.dirname(__file__)
    galaxia_params = test_filepath + '/galaxyModelParams_PopSyCLEv3.txt'

    synthetic.write_galaxia_params(output_root=output_root,
                                   longitude=1.25,
                                   latitude=-2.65,
                                   area=0.0001,
                                   seed=seed)

    synthetic.run_galaxia(output_root=output_root,
                          longitude=1.25,
                          latitude=-2.65,
                          area=0.0001,
                          galaxia_galaxy_model_filename=galaxia_params,
                          seed=seed)

    return output_root

@pytest.fixture()
def srun_popsyn(srun_galaxia):
    seed = 1

    ebf_file = srun_galaxia + '.ebf'
    output_root = srun_galaxia


    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='Raithel18',
                              overwrite=True,
                              seed=seed)

    return output_root

@pytest.fixture()
def srun_calc_events(srun_popsyn):
    seed = 1

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

@pytest.fixture()
def srun_refine_events(srun_calc_events):
    seed = 1

    input_root = srun_calc_events

    synthetic.refine_events(input_root=input_root,
                            filter_name='I',
                            photometric_system='ubv',
                            red_law='Damineli16',
                            overwrite=True,
                            output_file='default')

    output_root = input_root + '_refined_events_ubv_I_Damineli16'

    return output_root

@pytest.fixture
def mrun_galaxia():
    seed = 1

    output_root = 'data_test/test_Mrun'

    test_filepath = os.path.dirname(__file__)
    galaxia_params = test_filepath + '/galaxyModelParams_PopSyCLEv3.txt'

    synthetic.write_galaxia_params(output_root=output_root,
                                   longitude=1.25,
                                   latitude=-2.65,
                                   area=0.0001,
                                   seed=seed)

    synthetic.run_galaxia(output_root=output_root,
                          longitude=1.25,
                          latitude=-2.65,
                          area=0.0001,
                          galaxia_galaxy_model_filename=galaxia_params,
                          seed=seed)

    return output_root

@pytest.fixture()
def mrun_popsyn(mrun_galaxia):
    seed = 1

    ebf_file = mrun_galaxia + '.ebf'
    output_root = mrun_galaxia

    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)

    synthetic.perform_pop_syn(ebf_file=ebf_file,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='Raithel18',
                              multiplicity=multi_obj,
                              overwrite=True,
                              seed=seed)

    return output_root

@pytest.fixture()
def mrun_calc_events(mrun_popsyn):
    seed = 1

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

@pytest.fixture()
def mrun_refine_events(mrun_calc_events):
    input_root = mrun_calc_events

    synthetic.refine_events(input_root=input_root,
                            filter_name='I',
                            photometric_system='ubv',
                            red_law='Damineli16',
                            hdf5_file_comp=input_root + '_companions.h5',
                            overwrite=True,
                            output_file='default')

    output_root = input_root + '_refined_events_ubv_I_Damineli16'

    return output_root

@pytest.fixture()
def mrun_refine_binary(mrun_refine_events):
    input_root = mrun_refine_events
    events_prim_fits = input_root + '.fits'
    events_comp_fits = input_root + '_companions.fits'
    synthetic.refine_binary_events(events_prim_fits,
                                   events_comp_fits,
                                   filter_name='I',
                                   photometric_system='ubv', overwrite=True,
                                   output_file='default', save_phot=False)

    output_root = input_root + '_companions_rb'

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
def time_calc_event_time_loop():
    t1 = test_calc_event_time_loop(1)
    t4 = test_calc_event_time_loop(4)
    t8 = test_calc_event_time_loop(8)
    t12 = test_calc_event_time_loop(12)

    t_proc = np.array([t1, t4, t8, t12])
    t_thry = t_proc[0] * np.array([1, 1 / 4.0, 1 / 8.0, 1 / 12.0])
    n_proc = np.array([1, 4, 8, 12])

    fig = plt.figure(1, figsize=(6, 6))
    plt.clf()
    plt.plot(n_proc, t_thry, '.', label='Ideal')
    plt.plot(n_proc, t_proc, '.', label='Reality')
    plt.xlabel('N processors')
    plt.ylabel('Time (sec)')
    plt.legend()
    plt.show()

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


# FIXME
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
        
        # Makes mass slightly off from the star_masses
        mass = np.linspace(0.1, 5.1, 10)
        
        # All odd indices are Multiple systems
        isMulti = np.empty((len(mass),))
        isMulti[::2] = False
        isMulti[1::2] = True
        
        # Makes first two white dwarfs and the rest stars
        phase = np.empty((len(mass),))
        phase[0:2] = 101
        phase[2:] = 1.0
        
        star_systems = Table([mass, isMulti, phase],
                             names=['mass', 'isMultiple', 'phase'])
                             
        system_idx = [1,1,1,3,5,5,7,9]
        comp_phase = np.ones((len(system_idx)), dtype=int)
        comp_mass = np.random.rand(len(system_idx)) * mass[system_idx]  # uniform mass ratio.
        companions = Table([system_idx, comp_phase, comp_mass],
                             names=['system_idx', 'phase', 'mass'])
        
    cluster_check = FakeCluster()
        
    companion_check = synthetic._add_multiples(star_masses_check, cluster_check)
        
    # 1,1,1 should be eliminated because the primary is a WD
    # 9 should be eliminated because it's too big
    # 3, 5, 5, and 7 should be pointed to their star_masses counterpart
    system_idx_correct = [5,8,8,11]

    if not np.array_equal(system_idx_correct, companion_check['system_idx']):
        raise Exception("_add_multiples() is behaving unexpectedly")
    
    return

def test_system_mass(mrun_popsyn):
    """
    Test whether the system mass (current) is really the sum of the individual component masses.

    Run working_test_all_Mrun() first.
    """

    test_filepath = os.path.dirname(__file__)
    seed = 1

    galaxia_input = test_filepath + '/' + 'test.ebf'
    output_root = test_filepath + '/' + 'test_system_mass'

    synthetic.perform_pop_syn(ebf_file=galaxia_input,
                              output_root=output_root,
                              iso_dir='/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number=None,
                              BH_kick_speed_mean=100,
                              NS_kick_speed_mean=350,
                              IFMR='SukhboldN20',#'Raithel18',
                              overwrite=True,
                              multiplicity=multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2),
                              seed=seed)

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
        for ss in range(1000):
            primary = tab_syst[ss]

            if primary['isMultiple']:
                cdx = np.where(tab_comp['system_idx'] == primary['obj_id'])[0]
                companions = tab_comp[cdx]

                # Sum up the primary + companion current masses.
                mass = np.sum(companions['mass_current'])
                mass += primary['mass']

                # Sum up the primary + companion zams masses.
                zmass = np.sum(companions['zams_mass'])
                zmass += primary['zams_mass']

                # While we are here, check the current mass is < the zams mass.
                # Adding + 10**-4 since at very low masses sometimes isochrone
                # puts zams_mass < mass
                assert primary['mass'] <= primary['zams_mass'] + 10**-4
                for comp in companions:
                    assert comp['mass_current'] <= comp['zams_mass'] + 10**-4
                    
                #FIXME: add a check for companion < primary mass

                # Confirm that the sum of masses equals the system mass.
                assert isclose(mass, primary['systemMass'], abs_tol=1e-5)

                # Confirm that the current system mass is less than the zams system mass.
                assert mass < zmass

                # Keep track of how many we checked.
                n_comp_checked += 1

        print(f'Field {field}: Successfully checked system masses of {n_comp_checked} multiple systems.')
    return

# FIXME
def test_single_CO_frac():
    """
    Checks that the CO fraction of objects greater than 0.1 Msun
    is about 8.2%
    """
    test_filepath = os.path.dirname(__file__)
    seed = 1
    
    synthetic.perform_pop_syn(ebf_file = test_filepath + '/' + 'test.ebf',
                              output_root = test_filepath + '/' + 'test',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None,
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              overwrite=True,
                              seed = seed)
    test_hdf5 = h5py.File(test_filepath + '/' + 'test.h5', 'r')
    lower_mass_cutoff = 0.1 #Msun
    CO_frac = calc_CO_frac_mass_cutoff(test_hdf5, lower_mass_cutoff)
    
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
    return  CO_frac

# FIXME
def test_multiplicity_properties():
    """
    Checks that the multiplicity fraction of objects > 0.5 Msun is about 47%
    and that the minimum semimajor axis is greater than 10^-2
    """
    test_filepath = os.path.dirname(__file__)
    seed = 1
    
    synthetic.perform_pop_syn(ebf_file = test_filepath + '/' + 'test.ebf',
                              output_root = test_filepath + '/' + 'test_Mrun',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None,
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              overwrite=True,
                              multiplicity=multiplicity.MultiplicityResolvedDK(companion_max = True),
                              seed = seed)
    test_hdf5 = h5py.File(test_filepath + '/' + 'test_Mrun.h5', 'r')
    lower_mass_cutoff = 0.5 #Msun
    multiplicity_frac = calc_multiplicity_frac_mass_cutoff(test_hdf5, lower_mass_cutoff)
    
    
    precalc_mult_frac = 0.4727
    precalc_mult_number = 2052
    precalc_total_number = 4341
    precalc_error = precalc_mult_frac*np.sqrt((np.sqrt(precalc_mult_number)/precalc_mult_number)**2 + (np.sqrt(precalc_total_number)/precalc_total_number)**2)
    
    assert(np.abs(multiplicity_frac - precalc_mult_frac) < precalc_error)
    
    test_hdf5.close()
    
    
    test_hdf5_comp = h5py.File(test_filepath + '/' + 'test_Mrun_companions.h5', 'r')
    min_log_semimajor_axis = calc_min_semimajor_axis(test_hdf5_comp)
    
    min_log_semimajor_axis_imposed = -2
    
    assert(min_log_semimajor_axis_imposed < min_log_semimajor_axis)
    
    return

def calc_multiplicity_frac_mass_cutoff(hdf5_file, lower_mass_cutoff):
    subfield_list = list(hdf5_file.keys())[1:-2]
    multiples = 0
    total = 0
    for field in subfield_list:
        array = hdf5_file[field]
        multiples += len(np.where((array['isMultiple'] == 1) & (array['mass'] > lower_mass_cutoff))[0])
        total += len(np.where(array['mass'] > lower_mass_cutoff)[0])
        del array
    multiple_frac = multiples/total
    return  multiple_frac

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



# FIXME
def test_binary_angles():
    test_filepath = os.path.dirname(__file__)
    
    if os.path.exists(test_filepath + '/' + 'test_Mrun_events.fits') == False:
        synthetic.calc_events(hdf5_file = test_filepath + '/' + 'test_Mrun.h5',
                              output_root2 = test_filepath + '/' + 'test_Mrun',
                              radius_cut = 2,
                              obs_time = 1000,
                              n_obs = 11,
                              theta_frac = 2,
                              blend_rad = 0.65,
                              hdf5_file_comp = test_filepath + '/' + 'test_Mrun_companions.h5',
                              overwrite = True,
                              n_proc = 1)
        
    synthetic.refine_events(input_root = test_filepath + '/' + 'test_Mrun', 
                        filter_name = 'I',
                        photometric_system = 'ubv',
                        red_law = 'Damineli16', 
                        hdf5_file_comp = test_filepath + '/' + 'test_Mrun_companions.h5',
                        overwrite = True, 
                        output_file = 'default')
    
    
    test_companions_table = Table.read(test_filepath + '/' + 'test_Mrun' + '_refined_events_ubv_I_Damineli16_companions' + '.fits')
    
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

    

def working_test_all_Srun():
    """
    Testing an Srun (singles only)
    from beginning to end. This will generate all files and 
    check them against pre-run files.
    """
    
    seed = 1
    
    synthetic.write_galaxia_params(output_root = 'test',
                                   longitude = 1.25,
                                   latitude = -2.65,
                                   area = 0.0001,
                                   seed = seed)
    
    synthetic.run_galaxia(output_root = 'test',
                          longitude = 1.25,
                          latitude = -2.65,
                          area = 0.0001,
                          galaxia_galaxy_model_filename= '/g/lu/code/galaxia/docs/galaxyModelParams_PopSyCLEv3.txt',
                          seed = seed)
    
    #run_galaxia_result = filecmp.cmp('test_correct.ebf', 'test.ebf', shallow = False)
    #assert run_galaxia_result
    
    synthetic.perform_pop_syn(ebf_file = 'test.ebf',
                              output_root = 'test',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None,
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              overwrite=True,
                              seed = seed)
    
    perform_pop_syn_result = filecmp.cmp('test_correct.h5', 'test.h5', shallow = False) 
    #assert perform_pop_syn_result
    
    synthetic.calc_events(hdf5_file = 'test.h5',
                          output_root2 = 'test',
                          radius_cut = 2,
                          obs_time = 1000,
                          n_obs = 11,
                          theta_frac = 2,
                          blend_rad = 0.65,
                          overwrite = True,
                          n_proc = 1)
    
    calc_events_result = filecmp.cmp('test_correct.fits', 'test.fits', shallow = False) 
    #assert calc_events_result
    
    synthetic.refine_events(input_root = 'test',
                            filter_name = 'I',
                            photometric_system = 'ubv',
                            red_law = 'Damineli16',
                            overwrite = True,
                            output_file = 'default')
    
    refine_events_result = filecmp.cmp('test_correct_refined_events_ubv_I_Damineli16.fits', 'test_refined_events_ubv_I_Damineli16.fits', shallow = False) 
    #assert refine_events_result
    
    os.remove("test_galaxia_params.txt")
    os.remove("test.ebf")
    os.remove("test_galaxia.log")
    os.remove("test.h5")
    os.remove("test_perform_pop_syn.log")
    os.remove("test.fits")
    os.remove("test_blends.fits")
    os.remove("test_label.fits")
    os.remove("test_calc_events.log")
    os.remove("test_refined_events_ubv_I_Damineli16.fits")
    os.remove("test_refined_events_ubv_I_Damineli16.log")
    
    return

def working_test_all_Mrun():
    """
    Testing an Mrun (singles and multiples)
    from beginning to end. This will generate all files and 
    check them against pre-run files.

    Takes ~3 minutes to run.
    """
    
    seed = 1
    
    synthetic.write_galaxia_params(output_root = 'test_Mrun',
                                   longitude = 1.25,
                                   latitude = -2.65,
                                   area = 0.0001,
                                   seed = seed)
    
    synthetic.run_galaxia(output_root = 'test_Mrun',
                          longitude = 1.25,
                          latitude = -2.65,
                          area = 0.0001,
                          galaxia_galaxy_model_filename= '/g/lu/code/galaxia/docs/galaxyModelParams_PopSyCLEv3.txt',
                          seed = seed)
    
    run_galaxia_result = filecmp.cmp('test_correct_Mrun.ebf', 'test_Mrun.ebf', shallow = False) 
    assert run_galaxia_result 
    
    
    synthetic.perform_pop_syn(ebf_file = 'test_Mrun.ebf',
                              output_root = 'test_Mrun',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None,
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              multiplicity=multiplicity.MultiplicityResolvedDK(companion_max = True),
                              overwrite=True,
                              seed = seed)
    
    perform_pop_syn_result = filecmp.cmp('test_correct_Mrun.h5', 'test_Mrun.h5', shallow = False) 
    assert perform_pop_syn_result 
    
    perform_pop_syn_companions_result = filecmp.cmp('test_correct_Mrun_companions.h5', 'test_Mrun_companions.h5', shallow = False)
    assert perform_pop_syn_companions_result
    
    
    synthetic.calc_events(hdf5_file = 'test_Mrun.h5',
                          output_root2 = 'test_Mrun',
                          radius_cut = 2,
                          obs_time = 1000,
                          n_obs = 11,
                          theta_frac = 2,
                          blend_rad = 0.65,
                          hdf5_file_comp = 'test_Mrun_companions.h5',
                          overwrite = True,
                          n_proc = 1)
    
    calc_events_result = filecmp.cmp('test_correct_Mrun.fits', 'test_Mrun.fits', shallow = False) 
    assert calc_events_result 
    
    
    synthetic.refine_events(input_root = 'test_Mrun', 
                        filter_name = 'I',
                        photometric_system = 'ubv',
                        red_law = 'Damineli16', 
                        hdf5_file_comp = 'test_Mrun_companions.h5',
                        overwrite = True, 
                        output_file = 'default')
    
    refine_events_result = filecmp.cmp('test_correct_Mrun_refined_events_ubv_I_Damineli16.fits', 
                                       'test_Mrun_refined_events_ubv_I_Damineli16.fits', shallow = False) 
    assert refine_events_result 
    
    refine_events_companions_result = filecmp.cmp('test_correct_Mrun_refined_events_ubv_I_Damineli16_companions.fits',
                                              'test_Mrun_refined_events_ubv_I_Damineli16_companions.fits', shallow = False) 
    assert refine_events_companions_result 
    
    synthetic.refine_binary_events('test_Mrun_refined_events_ubv_I_Damineli16.fits', 
                                   'test_Mrun_refined_events_ubv_I_Damineli16_companions.fits', filter_name = 'I',
                                   photometric_system = 'ubv', overwrite = True,
                                   output_file = 'default', save_phot = False)
    
    refine_binary_events_result = filecmp.cmp('test_Mrun_correct_refined_events_ubv_I_Damineli16_companions_rb.fits',
                                          'test_Mrun_refined_events_ubv_I_Damineli16_companions_rb.fits', shallow = False) 
    assert refine_binary_events_result
    
    refine_binary_events_multi_peaks_result = filecmp.cmp('test_Mrun_correct_refined_events_ubv_I_Damineli16_companions_rb_multi_peaks.fits',
                                                      'test_Mrun_refined_events_ubv_I_Damineli16_companions_rb_multi_peaks.fits', shallow = False) 
    assert refine_binary_events_multi_peaks_result
    
    return
    
    
def generate_Srun_files():
    """
    This generates the correct Srun (singles only)
    files. This should !ONLY! be run if the outputs of the code
    were expected to change from when they were originally made!!!

    Takes ~3 minutes to run.
    """
    
    seed = 1
    
    synthetic.write_galaxia_params(output_root = 'test_correct',
                                   longitude = 1.25,
                                   latitude = -2.65,
                                   area = 0.0001,
                                   seed = seed)
    
    synthetic.run_galaxia(output_root = 'test_correct',
                          longitude = 1.25,
                          latitude = -2.65,
                          area = 0.0001,
                          galaxia_galaxy_model_filename= '/g/lu/code/galaxia/docs/galaxyModelParams_PopSyCLEv3.txt',
                          seed = seed)
    
    
    synthetic.perform_pop_syn(ebf_file = 'test_correct.ebf',
                              output_root = 'test_correct',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None, 
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              overwrite=False,
                              seed = seed)
    
    
    synthetic.calc_events(hdf5_file = 'test_correct.h5', 
                          output_root2 = 'test_correct', 
                          radius_cut = 2, 
                          obs_time = 1000, 
                          n_obs = 11, 
                          theta_frac = 2, 
                          blend_rad = 0.65, 
                          overwrite = False, 
                          n_proc = 1)
    
    
    synthetic.refine_events(input_root = 'test_correct', 
                            filter_name = 'I',
                            photometric_system = 'ubv',
                            red_law = 'Damineli16', 
                            overwrite = False, 
                            output_file = 'default')
    
    return
    
    
def generate_Mrun_files():
    """
    This generates the correct Mrun (singles and multiples)
    files. This should !ONLY! be run if the outputs of the code
    were expected to change from when they were originally made!!!

    Takes ~4 minutes to run.
    """
    
    seed = 1
    
    synthetic.write_galaxia_params(output_root = 'test_correct_Mrun',
                                   longitude = 1.25,
                                   latitude = -2.65,
                                   area = 0.0001,
                                   seed = seed)
    
    synthetic.run_galaxia(output_root = 'test_correct_Mrun',
                          longitude = 1.25,
                          latitude = -2.65,
                          area = 0.0001,
                          galaxia_galaxy_model_filename= '/g/lu/code/galaxia/docs/galaxyModelParams_PopSyCLEv3.txt',
                          seed = seed)
    
    
    synthetic.perform_pop_syn(ebf_file = 'test_correct_Mrun.ebf',
                              output_root = 'test_correct_Mrun',
                              iso_dir = '/g/lu/models/PopSyCLE_isochrones',
                              bin_edges_number = None,
                              BH_kick_speed_mean = 100,
                              NS_kick_speed_mean = 350,
                              IFMR = 'Raithel18',
                              multiplicity=multiplicity.MultiplicityResolvedDK(companion_max = True),
                              overwrite=False,
                              seed = seed)
    
    
    synthetic.calc_events(hdf5_file = 'test_correct_Mrun.h5',
                          output_root2 = 'test_correct_Mrun',
                          radius_cut = 2,
                          obs_time = 1000,
                          n_obs = 11,
                          theta_frac = 2,
                          blend_rad = 0.65,
                          hdf5_file_comp = 'test_correct_Mrun_companions.h5',
                          overwrite = False,
                          n_proc = 1)
    
    
    synthetic.refine_events(input_root = 'test_correct_Mrun', 
                        filter_name = 'I',
                        photometric_system = 'ubv',
                        red_law = 'Damineli16', 
                        hdf5_file_comp = 'test_correct_Mrun_companions.h5',
                        overwrite = False, 
                        output_file = 'default')
    
    synthetic.refine_binary_events('test_correct_Mrun_refined_events_ubv_I_Damineli16.fits', 
                                   'test_correct_Mrun_refined_events_ubv_I_Damineli16_companions.fits', filter_name = 'I',
                                   photometric_system = 'ubv', overwrite = False,
                                   output_file = 'default', save_phot = False)
    
    return
