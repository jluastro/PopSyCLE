.. _howto:

====================
Installation
====================
Fix this...

.. _Dependencies:

Dependencies
--------------
In addition to installing PopSyCLE, the user will need 

* Galaxia
* SPISEA (PyPopStar)
* Other python stuff...

.. _Installing_Galaxia:

Installing Galaxia
------------------
https://github.com/jluastro/galaxia

.. _Installing_SPISEA:

Installing SPISEA
-----------------
THIS NEEDS TO BE UPDATED https://github.com/astropy/PopStar

.. _Installing_Python_Libraries:

Installing Python Libraries
----------------------------
We recommend the Anaconda distribution.
In particular, numpy v1.13 or higher is required, along with Astropy and H5py.


===============
Reading files
===============
PopSyCLE use all sorts of different file formats. 
It can easily get confusing, so here is a short guide to the basics.

.. HDF5:

How to read HDF5 files
----------------------
Within the HDF5 file are datasets that store the information. 
It is kind of like a dictionary in python-- the dataset can be manipulated just like a numpy array.

CHECK: Sam/Michael is this how it works for HDF5 still?

First, go to the directory containing the HDF5 file you want to open. Next, start ipython. Then type the following::

    import h5py
    hf = h5py.File(`filename.h5', `r')

If you want to see the names of all the datasets in an HDF5 file, type the following::

   list(hf.keys())

Suppose you want to work with the dataset named ``dname``.
To access the dataset, type::

   dset = hf[`dname']

Note that only one person at a time can work on an open HDF5 file. Thus, at the end, you need to close the file::

     hf.close()

.. EBF:

How to read EBF files
---------------------
The EBF file is basically a dictionary in python. 
The output of Galaxia is in the EBF format. 

First, go to the directory containing the EBF file you want to open. 
Next, start ipython. 
Then type the following::

     import ebf 
     ef = ebf.read(`filename.ebf', `/')

If you want to see the names of all the keys in the EBF file, type the following::

   ef.keys()

Suppose you want to work with the key ``xkey``. To access that part of the file, type::

        x = ef[`xkey']

Now ``x`` is just a numpy array and can manipulated as such.

.. FITS:

How to read FITS table files
----------------------------
First, go to the directory containing the fits file you want to open. Next, start ipython, Then type the following::

       from astropy.table import Table
       tab = Table.read(`table.fits')

To view the entire table, just type ``tab``. The table works similar to a python dictionary. 
The column names are the keys of the dictionary, and the dictionary name in this case is ``tab``.

To view the header information/metadata, type ``tab.meta``.

\subsubsection{Stars and compact objects (HDF5)}

The data output contained in the HDF5 datasets are a combination of outputs that come directly from Galaxia, and outputs we ourselves have calculated or defined.

CHECK: PHOTOMETRIC SYSTEM KEYS... CAPTIAL OR LOWERCASE? (Michael?)

+----------+----------------------------------------------+--------------+
|Tag name  |Brief Description                             |Units         |
+==========+==============================================+==============+
|zams_mass |ZAMS mass                                     |Msun          |
+----------+----------------------------------------------+--------------+
|rem_id    |Integer indicating the remnant object type    |N/A           |
|          |(more details in tag description)             |              |
+----------+----------------------------------------------+--------------+
|obj_id    |Object ID-- unique integer to identify        |N/A           |
|          |star/compact object                           |              |
+----------+----------------------------------------------+--------------+
|mass      |Current mass                                  |Msun          |
+----------+----------------------------------------------+--------------+
|px        |Heliocentric x position                       |kpc           |
+----------+----------------------------------------------+--------------+
|py        |Heliocentric y position                       |kpc           | 
+----------+----------------------------------------------+--------------+
|pz        |Heliocentric z position                       |kpc           |
+----------+----------------------------------------------+--------------+
|vx        |Heliocentric x velocity                       |km/s          |
+----------+----------------------------------------------+--------------+
|vy        |Heliocentric y velocity                       |km/s          |
+----------+----------------------------------------------+--------------+
|vz        |Heliocentric z velocity                       |km/s          |
+----------+----------------------------------------------+--------------+
|rad       |Galactic radial distance                      |kpc           |
+----------+----------------------------------------------+--------------+
|glat      |Galactic latitude                             |deg           |
+----------+----------------------------------------------+--------------+
|glon      |Galactic longitude                            |deg           |
+----------+----------------------------------------------+--------------+
|vr        |Galactic radial velocity                      |km/s          |
+----------+----------------------------------------------+--------------+
|mu_b      |Galactic proper motion, b component           |mas/yr        |
+----------+----------------------------------------------+--------------+
|mu_lcosb  |Galactic proper motion, l component           |mas/yr        |
+----------+----------------------------------------------+--------------+
|age       |Age                                           |log(age/yr)   | 
+----------+----------------------------------------------+--------------+
|popid     |Population ID-- integer indicating the        |N/A           |
|          |population type ranging from 0 to 9           |              |
+----------+----------------------------------------------+--------------+
|ubv_x     |UBV photometric system,                       |mag           |
|          |x = {U, B, V, R, I, J, H, K}                  |              |
|          |absolute magnitude                            |              |
+----------+----------------------------------------------+--------------+
|exbv      |Extinction E(B-V) at the location of star     |mag           |
|          |given by 3-D Schlegel extinction maps         |              |
+----------+----------------------------------------------+--------------+
|teff      |Effective temperature                         |log(T/Kelvin) | 
+----------+----------------------------------------------+--------------+
|grav      |Surface gravity                               |log(gravity)  | 
+----------+----------------------------------------------+--------------+
|mbol      |Bolometric magnitude                          |log(L/Lsun)   | 
+----------+----------------------------------------------+--------------+
|feh       |Metallicity                                   |[Fe/H]        | 
+----------+----------------------------------------------+--------------+


\subsubsection{Event candidates table (Astropy FITS table)}

The event candidates table is very similar to the HDF5 file created by perform_pop_syn.
The main difference is there is a lot less output, so instead of writing it into an HDF5 file, we use an Astropy table.
Additionally, many of the tag names have a parenthetical (_L, _S) after them.
This is to indicate that there is one tag for the lens (L) and one for the source (S), since for a given event, 
you need to have both a lens and a source, and each of them has a mass, velocity, position, etc.
For example, zams_mass_L is the ZAMS mass of the lens, and age_S is the log(age/yr) of the source. 

+--------------------+----------------------------------------------+--------------+
|Tag name            |Brief Description                             |Units         |
+====================+==============================================+==============+
|zams_mass (_L, _S)  |ZAMS mass                                     |Msun          |
+--------------------+----------------------------------------------+--------------+
|rem_id (_L, _S)     |Integer indicating the remnant object type    |N/A           |
|                    |(more details in tag description)             |              |
+--------------------+----------------------------------------------+--------------+
|obj_id (_L, _S)     |Object ID-- unique integer to identify        |N/A           |
|                    |star/compact object                           |              |
+--------------------+----------------------------------------------+--------------+
|mass (_L, _S)       |Current mass                                  |Msun          |
+--------------------+----------------------------------------------+--------------+
|px (_L, _S)         |Heliocentric x position                       |kpc           |
+--------------------+----------------------------------------------+--------------+
|py (_L, _S)         |Heliocentric y position                       |kpc           | 
+--------------------+----------------------------------------------+--------------+
|pz (_L, _S)         |Heliocentric z position                       |kpc           |
+--------------------+----------------------------------------------+--------------+
|vx (_L, _S)         |Heliocentric x velocity                       |km/s          |
+--------------------+----------------------------------------------+--------------+
|vy (_L, _S)         |Heliocentric y velocity                       |km/s          |
+--------------------+----------------------------------------------+--------------+
|vz (_L, _S)         |Heliocentric z velocity                       |km/s          |
+--------------------+----------------------------------------------+--------------+
|rad (_L, _S)        |Galactic radial distance                      |kpc           |
+--------------------+----------------------------------------------+--------------+
|glat (_L, _S)       |Galactic latitude                             |deg           |
+--------------------+----------------------------------------------+--------------+
|glon (_L, _S)       |Galactic longitude                            |deg           |
+--------------------+----------------------------------------------+--------------+
|vr (_L, _S)         |Galactic radial velocity                      |km/s          |
+--------------------+----------------------------------------------+--------------+
|mu_b (_L, _S)       |Galactic proper motion, b component           |mas/yr        |
+--------------------+----------------------------------------------+--------------+
|mu_lcosb (_L, _S)   |Galactic proper motion, l component           |mas/yr        |
+--------------------+----------------------------------------------+--------------+
|age (_L, _S)        |Age                                           |log(age/yr)   | 
+--------------------+----------------------------------------------+--------------+
|popid (_L, _S)      |Population ID-- integer indicating the        |N/A           |
|                    |population type ranging from 0 to 9           |              |
+--------------------+----------------------------------------------+--------------+
|ubv_x (_L, _S)      |UBV photometric system,                       |mag           |
|                    |x = {U, B, V, R, I, J, H, K}                  |              |
|                    |absolute magnitude                            |              |
+--------------------+----------------------------------------------+--------------+
|exbv (_L, _S)       |Extinction E(B-V) at the location of star     |mag           |
|                    |given by 3-D Schlegel extinction maps         |              |
+--------------------+----------------------------------------------+--------------+
|teff (_L, _S)       |Effective temperature                         |log(T/Kelvin) | 
+--------------------+----------------------------------------------+--------------+
|grav (_L, _S)       |Surface gravity                               |log(gravity)  | 
+--------------------+----------------------------------------------+--------------+
|mbol (_L, _S)       |Bolometric magnitude                          |log(L/Lsun)   | 
+--------------------+----------------------------------------------+--------------+
|feh (_L, _S)        |Metallicity                                   |[Fe/H]        | 
+--------------------+----------------------------------------------+--------------+
|theta_E             |(Angular) Einstein radius                     |mas           |
+--------------------+----------------------------------------------+--------------+
|mu_rel              |Relative source-lens proper motion            |mas/yr        |
+--------------------+----------------------------------------------+--------------+
|u0                  |Minimum source-lens separation during the     |dim'less      |
|                    |survey, normalized by theta_E                 |              | 
+--------------------+----------------------------------------------+--------------+
|t0                  |Time at which the minimum source-lens         |days          |
|                    |separation occurs                             |              | 
+--------------------+----------------------------------------------+--------------+
