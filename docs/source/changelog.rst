Change Log
==========

=========================
v2.1.0 (2026-04-01)
=========================
* Inclusion of Rubin and Roman filters by creating interpolations between UBV filters and those respective filters
    * Author: Sage Remulla in `PR 66 <https://github.com/jluastro/PopSyCLE/pull/66>`_
* Support for orbital motion in lightcurves generated and refactor of ``refine_binary_events()`` to use lightcurve.py functions and be compatible with updated BAGLE
    * Authors: Alina Hussain and Natasha Abrams
* Inclusion of SDSS filter support via Galaxia
    * Author: Natasha Abrams
* Support for multiple filters to be analyzed at once in ``refine_events()`` by using a ``filt_dict``
    * *NOTE* in run.py this deprecates using ``photometric_system`` and ``filter`` in ``run.py``, but it is still allowed when running ``refine_events()`` directly
    * Author: Sage Remulla in `PR 68 <https://github.com/jluastro/PopSyCLE/pull/68>`_ and `PR 75 <https://github.com/MovingUniverseLab/PopSyCLE/pull/75>`_
* Changes and tutorial for compatibility with `synthpop <https://github.com/synthpop-galaxy/synthpop>`_ as a replacement for galaxia+perform_pop_syn in some cases
    * Author: Macy Huston in `PR 67 <https://github.com/jluastro/PopSyCLE/pull/67>`_
* Black hole splitting to allow user to choose binary black hole fraction
    * Fixed bugs in base function
    * Added capability to run this via run.py and slurm scripts
    * Author: Natasha Abrams
* Bug Fixes and Misc
    * Make Einstein radius more precise by using astropy constants
    * Added some asserts to verify companions are matched between primary and companion .h5 files (just for debugging, there was no issue)
    * Updated test_synthetic.py multiplicity fraction and CO fraction. The code that produced these didn't change, but it appears that the numpy seed generation changed causing the fractions to shift by 1-2%.


=========================
v2.0.2 (2024-12-24)
=========================

* Functions to manipulate lightcurves
    * Generate BAGLE models and parameters given a table of events
* HDF5 analysis utilities
    * Return a table with lists of star systems and their RA, Dec, z position, and system apparent magnitude. This is useful for making stellar density maps, computing microlens event occurrence rates, etc.
    * Count number of stars in hdf5 file
    * Makes a fraction of black holes single after perform_pop_syn. This is a temporary fix since there is no binary evolution and all the black holes end up in binaries.
* Bug Fixes
    * Issues involving masked columns causing masked blends and magnitudes instead of nans or correctly added ones
    * Column names of magnitudes of companion table

=========================
v2.0.1 (2024-03-12)
=========================

Small utility bug fixes

=========================
v2.0.0 (2024-01-29)
=========================

New Features:
-------------
* Multiplicity
    * See `Abrams et al. 2025 <https://ui.adsabs.harvard.edu/abs/2025ApJ...980..103A/abstract>`_
    * See `new example <https://github.com/jluastro/PopSyCLE/blob/main/docs/PopSyCLE_example_multiples.ipynb>`_ for details of how to call functions to add multiple systems. 
        
        * run_galaxia() is not modified, but all other function calls can be optionally changed to add multiple systems.
	* New function refine_binary_events() added to simulate binary lightcurves using BAGLE.

* New metallicity-dependent Initial-Final Mass Relations (IFMR)
    * See `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_ for details.
    * You are now required to input an IFMR when calling ``perform_pop_syn`` by specifying ``IFMR = 'SukhboldN20'`` (the other options are ``‘Raithel18’`` and ``‘Spera15’``). The previous default was ``‘Raithel18’``.
* Maxwellian Kicks
    * A more sophisticated method for assigning birth kick velocities to newly formed compact objects was also added to PopSyCLE, described in more detail in Section 2.4 of `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_. Instead of a single set of birth kicks of 100 km/s for black holes and 350 km/s for neutron stars being applied, a Maxwellian distribution with means at those values respectively were applied.
    * NS masses drawn from distribution rather than a single mass
    * See Section 2.2 of `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_  for details.
* perform_pop_syn() multiprocessing
* Updated test suite to use pytest

Method Changes:
---------------
* Dynamic cluster generation
    * Previously, the initial-current mass ratio of the cluster was used to find the necessary input-mass of SPISEA cluster to match the final primary stellar mass and the Galaxia stellar mass, as described in Section 3 of `Lam et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...889...31L/abstract>`_. Now, a dynamic cluster generation method is used to build up enough mass. SPISEA clusters are generated in chunks of the matching Galaxia stellar mass and combined until enough final stellar mass is generated. The associated companions and interspersed compact objects are kept.
* u0 is now assigned a sign: 
    * u0 > 0 when source is to the east of the lens.
    * u0 < 0 when source is to the west of the lens.

Bug Fixes:
----------
* Δm fix
    * The calculation of Δm was corrected to take into account the proper neighbors of the system. Previously, this led to a stochastic issue which would lead to too many or too few systems being cut which would likely average out over a number of fields.

Misc:
-----
* Various memory and speed up fixes.
* All mentions of PyPopStar were converted to SPISEA to reflect the software name changes.

===================
v1.1.0 (2020-03-02)
===================
* Filter Systems, Reformatted HDF5 Files, Pipeline Execution
   * This release contains the abilities to add Zwicky Transient Facility photometry to the PopSyCLE data pipeline, with a structure set up for future support of additional filter systems. 
   * Additionally, hdf5 files generated by generate_pop_syn now contain the metallicity, bolometric luminosity, surface gravity, and effective temperature of stars which may be used for population selections and stellar typing. The structure of these hdf5 files have been reformatted to use numpy compound datatypes, allowing for easy user access without using the col_idx objects in previous releases. 
   * This release enables the execution of the PopSyCLE pipeline with the command line tool run.py, as well as support to run this pipeline executor on a Slurm scheduler.

===================
v1.0.0 (2019-09-28)
===================
* Initial release (`Lam et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...889...31L/abstract>`_).
