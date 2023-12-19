changelog
========

=========================
v2.0.0 (Insert Date Here)
=========================

New Features:
-------------
* Multiplicity
    * See Abrams et al., in prep for details
    * See new example HERE for details of how to call functions to add multiple systems. 
        
        * run_galaxia() is not modified, but all other function calls can be optionally changed to add multiple systems.
* New metallicity-dependent Initial-Final Mass Relations (IFMR)
    * See `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_ for details.
    * You are now required to input an IFMR when calling ``perform_pop_syn`` by specifying ``IFMR = 'SukhboldN20'`` (the other options are ``‘Raithel18’`` and ``‘Spera15’``). The previous default was ``‘Raithel18’``.
* Maxwellian Kicks
    * A more sophisticated method for assigning birth kick velocities to newly formed compact objects was also added to PopSyCLE, described in more detail in Section 2.4 of `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_. Instead of a single set of birth kicks of 100 km/s for black holes and 350 km/s for neutron stars being applied, a Maxwellian distribution with means at those values respectively were applied.
    * NS masses drawn from distribution rather than a single mass
    * See Section 2.2 of `Rose et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..116R/abstract>`_  for details.
* perform_pop_syn multiprocessing
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

==================
v1.0.0 (1-22-2020)
==================
* Initial release (`Lam et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...889...31L/abstract>`_).
