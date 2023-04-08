from popsycle import synthetic
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.coordinates as coord
from astropy.table import Table
from astropy import units
import numpy as np
import h5py
from multiprocessing import Pool
import itertools
import matplotlib.pyplot as plt
import time
import ebf
import astropy.units as u
import astropy.coordinates as coord

masyr_to_degday = 1.0 * (1.0e-3 / 3600.0) * (1.0 / 365.25)
kms_to_kpcday = 1.0 * (3.086 * 10 ** 16) ** -1 * 86400.0
au_to_kpc = 4.848 * 10 ** -9


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


def test_calc_events():
    test_data_dir = '/u/casey/scratch/work/microlens/galaxia_test/OGLE672/'
    hdf5_file = test_data_dir + 'OGLE672.h5'
    output_root2 = 'OGLE672_TEST'
    obs_time = 1000
    n_obs = 101
    radius_cut = 2
    theta_frac = 2
    blend_rad = 0.8
    microlens_path = '/u/casey/scratch/code/PopSyCLE'
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

def _test_add_multiples():
    """
    Checks that _add_multiples() is making the companions point to the correct Galaxia objects
    and eliminates the correct elements.
    """
    # Fake version of Galaxia star masses
    star_masses_check = np.linspace(0, 5, 10).tolist()
    
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
        companions = Table([system_idx],
                             names=['system_idx'])
        
    cluster_check = FakeCluster()
    
    companion_check = _add_multiples(star_masses_check, cluster_check)
        
    # 1,1,1 should be eliminated because the primary is a WD
    # 9 should be eliminated because it's too big
    # 3, 5,5, and 7 should be pointed to their star_masses counterpart
    system_idx_correct = [4,6,6,8]
    correct_companion_check = Table([system_idx_correct],
                                   names=['system_idx'])
    
    if not np.array_equal(correct_companion_check, companion_check):
        raise Exception("_add_multiples() is behaving unexpectedly")
    
    return

