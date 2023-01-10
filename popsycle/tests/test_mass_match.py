import pdb

import numpy as np
import pylab as plt
import pytest
from spisea.imf import multiplicity

from popsycle import synthetic
import time


@pytest.fixture
def mass_arrays():
    masses = {}

    # Work in ranges of log mass to better approximate IMF.
    lo_s = -0.5
    lo_g = -0.95
    hi_s = 2.0
    hi_g = 2.0
    sz_s = 1000
    sz_g = 3000

    logm_s = np.random.uniform(low=lo_s, high=hi_s, size=sz_s)
    logm_g = np.random.uniform(low=lo_g, high=hi_g, size=sz_g)

    masses['spisea'] = 10 ** logm_s
    masses['galaxia'] = 10 ** logm_g

    return masses

def make_mass_arrays(n_spisea=100, n_galaxia=300):
    masses = {}

    # Work in ranges of log mass to better approximate IMF.
    lo_s = -0.5
    lo_g = -0.95
    hi_s = 2.0
    hi_g = 1.9
    sz_s = n_spisea
    sz_g = n_galaxia

    logm_s = np.random.uniform(low=lo_s, high=hi_s, size=sz_s)
    logm_g = np.random.uniform(low=lo_g, high=hi_g, size=sz_g)

    masses['spisea'] = 10 ** logm_s
    masses['galaxia'] = 10 ** logm_g

    return masses

def make_mass_arrays_edge_cases():
    masses = {}

    masses['spisea'] = np.array([10.0, 0.7, 1.2, 2.0, 2.0, 2.7])
    masses['galaxia'] = np.array([9.0, 0.6, 1.4, 1.5, 2.5, 2.5, 3.0])

    return masses

def test_match_companions_kdtree(mass_arrays):
    m_s = mass_arrays['spisea']
    m_g = mass_arrays['galaxia']

    # Pre-filter any SPISEA stars that are more massive than the most-massive galaxia star.
    # We do this in synthetic in earlier steps. Replicate here for timing accuracy.
    gdx = np.where(m_s <= m_g.max())[0]
    m_s = m_s[gdx]

    t0 = time.time()

    # Run the traditional mass matching.
    closest_index_arr, mass_diff = synthetic._match_companions_kdtree(m_g, m_s)
    dt = time.time() - t0
    print()
    print(f'Runtime = {time.time() - t0:.3} sec')
    print(f'Number of matches = {len(closest_index_arr)} of {len(m_s)}')

    # count the number of duplicates that remain
    _, counts = np.unique(closest_index_arr, return_counts=True)
    n_duplicates = np.sum(counts > 1)
    print(f'Number of duplicates = {n_duplicates}')

    return dt

def test_match_companions_kdtree_nonneg(mass_arrays, plot=False):
    m_s = mass_arrays['spisea']
    m_g = mass_arrays['galaxia']

    # Pre-filter any SPISEA stars that are more massive than the most-massive galaxia star.
    # We do this in synthetic in earlier steps. Replicate here for timing accuracy.
    gdx = np.where(m_s <= m_g.max())[0]
    m_s = m_s[gdx]

    t0 = time.time()
    closest_index_arr, closest_mass_diff = synthetic._match_companions_kdtree_nonneg(m_g, m_s, verbose=1)
    dt = time.time() - t0
    print()
    print(f'Runtime = {dt:.3} sec')

    ####
    # Results
    ####
    good_matches = np.where(closest_index_arr != -1)[0]
    bad_matches = np.where(closest_index_arr == -1)[0]
    m_g_matched = np.ones_like(m_s)
    m_g_matched[good_matches] = m_g[closest_index_arr[good_matches]]
    m_g_matched[bad_matches] = None
    print(f'Number of matches = {len(good_matches)} of {len(m_s)}')

    assert len(bad_matches) < len(good_matches)
    np.testing.assert_array_less(m_s[good_matches], m_g_matched[good_matches])
    np.testing.assert_array_less(-0.0001, closest_mass_diff[good_matches])
    # np.testing.assert_almost_equal(closest_mass_diff[good_matches],
    #                                m_g_matched[good_matches] - m_s[good_matches],
    #                                decimal=2)

    # count the number of duplicates that remain
    foo, counts = np.unique(closest_index_arr[good_matches], return_counts=True)
    n_duplicates = np.sum(counts > 1)
    print(f'Number of duplicates = {n_duplicates}', foo[counts > 1])
    # pdb.set_trace()

    if plot:
        import pylab as plt
        plt.figure()
        plt.plot(m_s, closest_mass_diff / m_s, 'k.')
        plt.xlabel('m_SPISEA')
        plt.ylabel('Fractional mass diff (Galaxia - SPISEA)')

    return dt

def test_match_companions_diff_array(mass_arrays):
    m_s = mass_arrays['spisea']
    m_g = mass_arrays['galaxia']

    # Pre-filter any SPISEA stars that are more massive than the most-massive galaxia star.
    # We do this in synthetic in earlier steps. Replicate here for timing accuracy.
    gdx = np.where(m_s <= m_g.max())[0]
    m_s = m_s[gdx]

    t0 = time.time()
    closest_index_arr, closest_mass_diff = synthetic._match_companions_diff_array(m_g, m_s)
    print(f'\n Runtime = {time.time() - t0:.3} sec')

    ####
    # Results
    ####
    good_matches = np.where(closest_index_arr != -1)[0]
    bad_matches = np.where(closest_index_arr == -1)[0]
    m_g_matched = np.ones_like(m_s)
    m_g_matched[good_matches] = m_g[closest_index_arr[good_matches]]
    m_g_matched[bad_matches] = None

    assert len(bad_matches) < len(good_matches)
    np.testing.assert_array_less(m_s[good_matches], m_g_matched[good_matches])
    np.testing.assert_array_less(0, closest_mass_diff[good_matches])
    np.testing.assert_almost_equal(closest_mass_diff[good_matches],
                                   m_g_matched[good_matches] - m_s[good_matches],
                                   decimal=3)

    # count the number of duplicates that remain
    _, counts = np.unique(closest_index_arr, return_counts=True)
    n_duplicates = np.sum(counts > 1)
    print(f'Number of matches = {len(good_matches)} of {len(m_s)}')
    print(f'Number of duplicates = {n_duplicates}')

    return

def test_add_multiples_match_with_kdtree():
    ####
    # Make fake data to work on.
    ####
    # Artificial cluster parameters.
    log_age = 9.0
    iso_dir = '/g/lu/models/PopSyCLE_isochrones'
    cluster_mass = 50000.0
    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)
    feh = 0.0
    seed = 1

    # Make a Galaxia cluster (just use SPISEA with no companions).
    # We will only use the mass column for now.
    # SPISEA cluster.
    foo = synthetic._make_cluster(iso_dir=iso_dir,
                                  log_age=log_age,
                                  currentClusterMass=cluster_mass,
                                  multiplicity=None,
                                  IFMR='SukhboldN20',
                                  feh=feh, seed=seed+10)
    cluster_g, _, _ = foo
    idx = np.where(cluster_g.star_systems['phase'] < 100)[0]
    m_curr_g = cluster_g.star_systems[idx]['mass'].tolist()
    m_zams_g = cluster_g.star_systems[idx]['zams_mass'].tolist()

    # masses = make_mass_arrays(n_galaxia=300)
    # m_g = masses['galaxia']

    foo2 = synthetic._make_cluster(iso_dir=iso_dir,
                                  log_age=log_age,
                                  currentClusterMass=np.sum(m_curr_g),
                                  multiplicity=multi_obj,
                                  IFMR='SukhboldN20',
                                  feh=feh, seed=seed)
    cluster_s, _, _ = foo2

    print(f'GALAXIA: Current cluster mass = {np.sum(m_curr_g):.1f} Msun')
    print(f'GALAXIA: Total number of primaries = {len(m_curr_g)}')
    print(f'SPISEA:  Current cluster mass = {cluster_s.star_systems["mass"].sum():.1f} Msun')
    print(f'SPISEA:  Total number of primaries = {len(cluster_s.star_systems)}')
    print(f'SPISEA:  Total number of companions = {len(cluster_s.companions)}')



    ###
    # Test _add_multiples
    ###
    t0 = time.time()
    modified_companions = synthetic._add_multiples_all(m_zams_g, cluster_s)
    print(f'Runtime = {time.time() - t0:.3} sec')

    assert len(modified_companions) > 0.5 * len(cluster_s.companions)

    # Check the number of companions in the galaxia companions table
    # is less than or equal to the number of spisea companions. Can't create
    # extra companions.
    assert len(modified_companions) < len(cluster_s.companions)

    # Print out the maximum mass difference between the primaries.
    # zams_mass_match_diff = galaxia mass - spisea mass.
    print(f'min max diff = {modified_companions["zams_mass_match_diff"].min():.2f}')
    print(f'max max diff = {modified_companions["zams_mass_match_diff"].max():.2f}')

    # mass diff should always be positive ideally. This highlights the problem...
    # lots of negative mass differences.
    pdx = np.where(modified_companions['zams_mass_match_diff'] >= 0)[0]
    ndx = np.where(modified_companions['zams_mass_match_diff'] < 0)[0]
    print(f'N pos diff = {len(pdx)}')
    print(f'N neg diff = {len(ndx)}')

    return

def test_add_multiples_match_with_kdtree_nonneg():
    ####
    # Make fake data to work on.
    ####
    # Artificial cluster parameters.
    log_age = 9.0
    iso_dir = '/g/lu/models/PopSyCLE_isochrones'
    cluster_mass = 50000.0
    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)
    feh = 0.0
    seed = 1

    # Make a Galaxia cluster (just use SPISEA with no companions).
    # We will only use the mass column for now.
    # SPISEA cluster.
    foo = synthetic._make_cluster(iso_dir=iso_dir,
                                  log_age=log_age,
                                  currentClusterMass=cluster_mass,
                                  multiplicity=None,
                                  IFMR='SukhboldN20',
                                  feh=feh, seed=seed + 10)
    cluster_g, _, _ = foo
    idx = np.where(cluster_g.star_systems['phase'] < 100)[0]
    m_curr_g = cluster_g.star_systems[idx]['mass'].tolist()
    m_zams_g = cluster_g.star_systems[idx]['zams_mass'].tolist()

    # masses = make_mass_arrays(n_galaxia=300)
    # m_g = masses['galaxia']

    foo2 = synthetic._make_cluster(iso_dir=iso_dir,
                                   log_age=log_age,
                                   currentClusterMass=np.sum(m_curr_g),
                                   multiplicity=multi_obj,
                                   IFMR='SukhboldN20',
                                   feh=feh, seed=seed)
    cluster_s, _, _ = foo2

    print(f'GALAXIA: Current cluster mass = {np.sum(m_curr_g):.1f} Msun')
    print(f'GALAXIA: Total number of primaries = {len(m_curr_g)}')
    print(f'SPISEA:  Current cluster mass = {cluster_s.star_systems["mass"].sum():.1f} Msun')
    print(f'SPISEA:  Total number of primaries = {len(cluster_s.star_systems)}')
    print(f'SPISEA:  Total number of companions = {len(cluster_s.companions)}')

    ###
    # Test _add_multiples
    ###
    t0 = time.time()
    modified_companions = synthetic._add_multiples_some_unmatched(m_zams_g, cluster_s, verbose=1)
    print(f'Runtime = {time.time() - t0:.3} sec')

    assert len(modified_companions) > 0.5 * len(cluster_s.companions)

    # Check the number of companions in the galaxia companions table
    # is less than or equal to the number of spisea companions. Can't create
    # extra companions.
    assert len(modified_companions) < len(cluster_s.companions)

    # Check that we have some duplicates (i.e. multiple companions matched
    # to a single primary).
    _, dup_count = np.unique(modified_companions['system_idx'], return_counts=True)
    assert dup_count.max() > 1

    # Make sure that the pointers back to Galaxia stars all make sense.
    for cc in range(len(modified_companions)):
        comp = modified_companions[cc]

        # Assert that the system idx exists and is positive.
        assert comp['system_idx'] >= 0

        # Primary star system from Galaxia
        prim_g = cluster_g.star_systems[idx][comp['system_idx']]

        # Check that companion mass ratio is q < 1 (zams mass).
        assert prim_g['zams_mass'] > comp['zams_mass']

    # Check the maximum mass difference between the primaries
    # isn't too large.
    # zams_mass_match_diff = galaxia mass - spisea mass.
    assert modified_companions["zams_mass_match_diff"].min() >= 0
    assert modified_companions["zams_mass_match_diff"].max() < 1.0 # msun

    # mass diff should always be positive ideally. This highlights the problem...
    # lots of negative mass differences.
    pdx = np.where(modified_companions['zams_mass_match_diff'] >= 0)[0]
    ndx = np.where(modified_companions['zams_mass_match_diff'] < 0)[0]

    assert len(pdx) > len(ndx)
    assert len(ndx) == 0

    return

def test_add_multiples_match_with_diff_array():
    ####
    # Make fake data to work on.
    ####
    # Artificial cluster parameters.
    log_age = 9.0
    iso_dir = '/g/lu/models/PopSyCLE_isochrones'
    cluster_mass = 50000.0
    multi_obj = multiplicity.MultiplicityResolvedDK(companion_max=True, CSF_max=2)
    feh = 0.0
    seed = 1

    # Make a Galaxia cluster (just use SPISEA with no companions).
    # We will only use the mass column for now.
    # SPISEA cluster.
    foo = synthetic._make_cluster(iso_dir=iso_dir,
                                  log_age=log_age,
                                  currentClusterMass=cluster_mass,
                                  multiplicity=None,
                                  IFMR='SukhboldN20',
                                  feh=feh, seed=seed + 10)
    cluster_g, _, _ = foo
    idx = np.where(cluster_g.star_systems['phase'] < 100)[0]
    m_curr_g = cluster_g.star_systems[idx]['mass'].tolist()
    m_zams_g = cluster_g.star_systems[idx]['zams_mass'].tolist()

    # masses = make_mass_arrays(n_galaxia=300)
    # m_g = masses['galaxia']

    foo2 = synthetic._make_cluster(iso_dir=iso_dir,
                                   log_age=log_age,
                                   currentClusterMass=np.sum(m_curr_g),
                                   multiplicity=multi_obj,
                                   IFMR='SukhboldN20',
                                   feh=feh, seed=seed)
    cluster_s, _, _ = foo2

    print(f'GALAXIA: Current cluster mass = {np.sum(m_curr_g):.1f} Msun')
    print(f'GALAXIA: Total number of primaries = {len(m_curr_g)}')
    print(f'SPISEA:  Current cluster mass = {cluster_s.star_systems["mass"].sum():.1f} Msun')
    print(f'SPISEA:  Total number of primaries = {len(cluster_s.star_systems)}')
    print(f'SPISEA:  Total number of companions = {len(cluster_s.companions)}')

    ###
    # Test _add_multiples
    ###
    t0 = time.time()
    modified_companions = synthetic._add_multiples_new(m_zams_g, cluster_s, verbose=0)
    print(f'Runtime = {time.time() - t0:.3} sec')

    assert len(modified_companions) > 0.5 * len(cluster_s.companions)

    # Check the number of companions in the galaxia companions table
    # is less than or equal to the number of spisea companions. Can't create
    # extra companions.
    assert len(modified_companions) < len(cluster_s.companions)

    # Check that we have some duplicates (i.e. multiple companions matched
    # to a single primary).
    _, dup_count = np.unique(modified_companions['system_idx'], return_counts=True)
    assert dup_count.max() > 1

    # Make sure that the pointers back to Galaxia stars all make sense.
    for cc in range(len(modified_companions)):
        comp = modified_companions[cc]

        # Assert that the system idx exists and is positive.
        assert comp['system_idx'] >= 0

        # Primary star system from Galaxia
        prim_g = cluster_g.star_systems[idx][comp['system_idx']]

        # Check that companion mass ratio is q < 1 (zams mass).
        assert prim_g['zams_mass'] > comp['zams_mass']

    # Check the maximum mass difference between the primaries
    # isn't too large.
    # zams_mass_match_diff = galaxia mass - spisea mass.
    assert modified_companions["zams_mass_match_diff"].min() >= 0
    assert modified_companions["zams_mass_match_diff"].max() < 1.0 # msun

    # mass diff should always be positive ideally. This highlights the problem...
    # lots of negative mass differences.
    pdx = np.where(modified_companions['zams_mass_match_diff'] >= 0)[0]
    ndx = np.where(modified_companions['zams_mass_match_diff'] < 0)[0]

    assert len(pdx) > len(ndx)
    assert len(ndx) == 0

    return

def time_vs_nstar_old(n_stars=None):
    if n_stars is None:
        n_stars = np.array([50, 100, 500, 1000, 5000, 7500, 10000, 25000, 50000, 75000, 100000])

    dt_all = np.zeros(len(n_stars), dtype=float)
    for nn in range(len(n_stars)):
        n_star = n_stars[nn]
        print('\n\n')
        print(f'** nn = {n_star}')
        mass_arrays = make_mass_arrays(n_spisea=int(n_star), n_galaxia=int(2*n_star))
        dt_all[nn] = test_match_companions_kdtree(mass_arrays)

    plt.figure(1)
    plt.clf()
    plt.loglog(n_stars, dt_all, 'k-X')
    plt.ylim(1e-3, 1e3)
    plt.xlabel('N stars')
    plt.ylabel('Runtime (s)')
    plt.title('old kdtree')
    plt.savefig('mass_match_kdtree_time_vs_nstar_old.png')

    return n_stars, dt_all

def time_vs_nstar(n_stars=None):
    if n_stars is None:
        n_stars = np.array([50, 100, 500, 1000, 5000, 7500, 10000, 25000, 50000, 75000, 100000])

    dt_all = np.zeros(len(n_stars), dtype=float)
    for nn in range(len(n_stars)):
        n_star = n_stars[nn]
        print('\n\n')
        print(f'** nn = {n_star}')
        mass_arrays = make_mass_arrays(n_spisea=int(n_star), n_galaxia=int(2*n_star))
        dt_all[nn] = test_match_companions_kdtree_nonneg(mass_arrays)

    plt.figure(2)
    plt.clf()
    plt.loglog(n_stars, dt_all, 'k-X')
    plt.ylim(1e-3, 1e3)
    plt.xlabel('N stars')
    plt.ylabel('Runtime (s)')
    plt.title('new kdtree')
    plt.savefig('mass_match_kdtree_time_vs_nstar_new.png')

    return n_stars, dt_all

def time_vs_nstar_compare():
    # n_stars = np.array([50, 100, 500, 1000, 5000, 7500, 10000, 25000])#, 50000, 75000, 100000])
    n_stars = np.array([50, 100, 500, 1000, 5000, 7500, 10000, 25000, 50000, 75000, 100000])

    n_stars, dt_old = time_vs_nstar_old(n_stars)
    n_stars, dt_new = time_vs_nstar(n_stars)

    dt_ratio = dt_new / dt_old

    plt.close(3)

    plt.figure(3, figsize=(6, 6))
    plt.subplots_adjust(top=0.98, right=0.98)

    ax1 = plt.subplot(2, 1, 1)
    plt.loglog(n_stars, dt_ratio, 'k-X')
    plt.xlabel('N stars')
    plt.ylabel('New / Old Runtime')

    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    plt.loglog(n_stars, dt_old, 'k-X', label='old')
    plt.loglog(n_stars, dt_new, 'r-X', label='new')
    plt.ylim(1e-3, 1e3)
    plt.xlim(2e1, 2e5)
    plt.legend()

    # 1-to-1 line scaled from first point.
    xlo = 2e1
    xhi = 2e5
    ylo = dt_new[0] * (xlo / n_stars[0])
    yhi = dt_new[0] * (xhi / xlo)
    print(xlo, xhi)
    print(ylo, yhi)
    plt.plot([xlo, xhi], [ylo, yhi], ls='--', marker=None, color='grey')

    plt.xlabel('N stars')
    plt.ylabel('Runtime (s)')

    plt.savefig('mass_match_kdtree_time_vs_nstar_compare.png')

    return
