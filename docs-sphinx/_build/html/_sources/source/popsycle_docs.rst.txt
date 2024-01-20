popsycle_docs
=============

.. contents:: Table of Contents
   :backlinks: none


==============
1 Installation
==============

1.1 Installing Galaxia
-----------------------

1.1.1 Installation
++++++++++++++++++

        Note there are also instructions from the creators of Galaxia at `<http://galaxia.sourceforge.net/Galaxia3pub.html>`_; 
        what follows is a very explicit version of that.
        
            #. Go to `<https://sourceforge.net/projects/galaxia/files/>`_ and download Galaxia by clicking the big green button.
            #. Go to your Downloads folder and untar the file by double-clicking it.
            #. Move the untar'd folder (which should be called something like galaxia-0.7.2) to your home directory. 
               (That's the directory where you get sent if you ``cd`` and don't put a location).
            #. In your home directory, make a directory called GalaxiaData (i.e. ``mkdir GalaxiaData``).
            #. Move to the galaxia-0.7.2 directory (i.e. ``cd galaxia-0.7.2``), and in there, do the following 
               (replacing ``/your/home/directory/`` with your home directory as appropriate):
                .. code-block:: bash

                    ./configure --datadir=/your/home/directory/GalaxiaData/
                    make
                    sudo make install
                    cp -r GalaxiaData/ /your/home/directory/GalaxiaData/
    
            #. Move back into your home directory (i.e. ``cd``), then run the following:
                .. code-block:: bash

                    galaxia -s warp
    
        That should be it! 
        You don't have to install ebf if you just download PopSyCLE!
        
        NOTE: the instructions in step 5 are assuming a root install. 
        If you want to do a local install, you need to have a folder for the software to be installed in.
        For example, in my home directory I made a `sw` directory (/your/home/directory/sw) for Galaxia to be installed in.
        You can then run the following instead:
            .. code-block:: bash
                
                ./configure --prefix=/your/home/directory/sw --datadir=/your/home/directory/GalaxiaData/
                make
                make install
                cp -r GalaxiaData/ /your/home/directory/GalaxiaData/

        You also need to export Galaxia to your path. 
        In your .bash_profile or .zshenv add the line ``export PATH=\$PATH:/your/home/directory/sw/bin``.
        Then proceed with step 6 in the installation instructions.

1.1.2 Uninstallation
++++++++++++++++++++

        You need to remove the compiled galaxia code (you can find where it is by typing ``which galaxia`` in the terminal), 
        the GalaxiaData directory, and you might as well remove the galaxia-0.7.2 directory also.
        When you do ``which galaxia`` nothing should be returned.

1.1.3 Parameter modification
++++++++++++++++++++++++++++

        Suppose you want to change the pattern speed in Galaxia.
        To do this, follow the installation instructions up to and including step 4.
        Then do the following:
            #. Move to the galaxia-0.7.2/src directory.
            #. Open the Population.h file with your favorite text editor.
            #. Find the pattern speed (in this case by searching for 71.62) and replace with your desired value (in this case 40.00).
            #. Save the change.
        Now return to step 5 in the installation instruction and proceed as instructed.

1.2 Installing SPISEA
----------------------

    SPISEA can be installed by cloning the repository from `<https://github.com/astropy/SPISEA>`_ and following the instructions.

1.3 Installing BAGLE
----------------------

    BAGLE can be installed by cloning the repository from `<https://github.com/MovingUniverseLab/BAGLE_Microlensing>`_ and following the instructions.

1.4 Installing Python libraries
--------------------------------

    We recommend the Anaconda distribution.
    In particular, numpy v1.24 or higher is required, along with Astropy and H5py.

===============
2 Reading Files
===============
PopSyCLE uses all sorts of different file formats. It can easily get
confusing, so here is a short guide to the basics.

2.1 How to read HDF5 files
---------------------------

    Within the HDF5 file are datasets that store the information. It is kind
    of like a dictionary in python-- the dataset can be manipulated just
    like a numpy array.
    
    First, go to the directory containing the HDF5 file you want to open.
    Next, start ipython. Then type the following:

    .. code-block:: python
    
        import h5py
        hf = h5py.File('filename.h5', 'r') 
    
    If you want to see the names of all the datasets in an HDF5 file, type
    the following:
    
    .. code-block:: python
    
        list(hf.keys())
    
    Suppose you want to work with the dataset named dname
    
       To access the dataset, type:

        .. code-block:: python
    
            dset = hf['dname']
    
    ..
    
       To view the columns of the dset, you can type:
    
        .. code-block:: python    
            
            dset.dtype.names
    
    ..
    
       To access a column in the database, you can use the column names
       listed below. i.e. for mass of all the objects in a dataset you can
       use:

        .. code-block:: python 
    
            dset_masses = hf['dname']['mass']
    
    ..
    
       Note that only one person at a time can work on an open HDF5 file.
       Thus, at the end, you need to close the file:
    
        .. code-block:: python 
            
            hf.close()
    
    ..

2.2 How to read EBF files
--------------------------
    The EBF file is basically a dictionary in python. The output of
    Galaxia is in the EBF format.
    
    First, go to the directory containing the EBF file you want to open.
    Next, start ipython. Then type the following:
    
    .. code-block:: python
        
        from popsycle import ebf
        ef = ebf.read('filename.ebf', '/')
    
    ..
    
    If you want to see the names of all the keys in the EBF file, type
    the following:
    
    .. code-block:: python
       
        ef.keys()
    
    ..
    
    Suppose you want to work with the key xkey . To access that part of
    the file, type:
    
    .. code-block:: python    

        x = ef['xkey']
    
    ..
    
    Now x is just a numpy array and can manipulated as such
    
    You can also access just that key from the beginning instead of
    loading in the entire ebf file by:
    
    .. code-block:: python     
        
        ef = ebf.read('filename.ebf', '/xkey')
    
    ..

2.3 How to read FITS table files
---------------------------------
    First, go to the directory containing the fits file you want to
    open. Next, start ipython, Then type the following:
    
    .. code-block:: python

        from astropy.table import Table
        tab = Table.read('table.fits')
    
    ..
    
    To view the entire table, just type tab . The table works similar to
    a python dictionary or like a pandas dataframe. The column names are
    the keys of the dictionary, and the dictionary name in this case is
    tab .
    
    To view the header information/metadata, type
    
    .. code-block:: python
        
        tab.meta
    
    ..
    
    To view the column names type
    
    .. code-block:: python
        
        tab.columns
    
    ..

==========================
3 Description of Pipeline
==========================

3.1 Pipeline without multiple systems
-------------------------------------

   To run without companions, first, run Galaxia to create an EBF file,
   which produces a synthetic survey, i.e. a bunch of stars. Next, run
   population synthesis (perform_pop_syn) to inject compact objects into
   the synthetic survey; both the compact objects and stars are saved in
   an HDF5 file. Then run a synthetic survey (calc_events and
   refine_events) that will produce a list of microlensing events, which
   are listed in a FITS file.

.. image:: popsycle_docs_images/media/pipeline.png
   :width: 3.2375in
   :height: 3.26528in
   :align: center

3.2 Pipeline with multiple systems
----------------------------------

   To run with companions (with changed steps marked in **bold**),
   first, run Galaxia to create an EBF file, which produces a synthetic
   survey, i.e. a bunch of stars. Next, run population synthesis
   (perform_pop_syn) to inject compact objects **and companions** into
   the synthetic survey. Both the compact objects and stars are saved in
   an HDF5 file **and the companions are stored in a separate hdf5
   file**. Then run a synthetic survey (calc_events and refine_events)
   that will produce a list of microlensing events, which are listed in
   a FITS file **with a separate FITS file for associated companions
   after refine_events**. **Then model the binary lens events
   (refine_binary_events) which will produce some additional
   characteristics from the lightcurves and a description of all the
   peaks, which are listed in two FITS files.**

.. image:: popsycle_docs_images/media/pipeline_w_multiples.png
   :align: center

==========
4 Outputs
==========

In addition to the outputs described below, each function produces
a text log file that lists the input parameters

4.1 perform_pop_syn
---------------------

4.1.1 Label/summary file (Astropy FITS Table)
+++++++++++++++++++++++++++++++++++++++++++++
        For now see

4.1.2 Stars and Compact Objects (primaries) (HDF5)
++++++++++++++++++++++++++++++++++++++++++++++++++
        The data output contained in the HDF5 datasets are a combination of
        outputs that come directly from Galaxia, and outputs we ourselves
        have calculated or defined.

       Default name: *root*.h5

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    zams_mass          |    ZAMS mass          |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    mass               |    Current mass       |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    systemMass         |    Sum of mass of     |    M⊙                 |
|                       |    primary and        |                       |
|                       |    companions (if     |                       |
|                       |    existent)          |                       |
+-----------------------+-----------------------+-----------------------+
|    px                 |    Heliocentric x     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    py                 |    Heliocentric y     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    pz                 |    Heliocentric z     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    vx                 |    Heliocentric x     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vy                 |    Heliocentric y     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vz                 |    Heliocentric z     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    age                |    Age                |    log(age/yr)        |
+-----------------------+-----------------------+-----------------------+
|    popid              |    Population ID -    |    N/A                |
|                       |    integer indicating |                       |
|                       |    the population     |                       |
|                       |    type ranging from  |                       |
|                       |    0 to 9 (see        |                       |
|                       |    Additional         |                       |
|                       |    Descriptions       |                       |
|                       |    below)             |                       |
+-----------------------+-----------------------+-----------------------+
|    exbv               |    Extinction E(B-V)  |    mag                |
|                       |    at the location of |                       |
|                       |    star given by 3-D  |                       |
|                       |    Schlegel           |                       |
|                       |    extinction maps    |                       |
+-----------------------+-----------------------+-----------------------+
|    glat               |    Galactic latitude  |    deg                |
+-----------------------+-----------------------+-----------------------+
|    glon               |    Galactic longitude |    deg                |
+-----------------------+-----------------------+-----------------------+
|    mbol               |    Bolometric         |    log(L/L⊙)          |
|                       |    magnitude          |                       |
+-----------------------+-----------------------+-----------------------+
|    grav               |    Surface gravity    |    log(gravity)       |
+-----------------------+-----------------------+-----------------------+
|    teff               |    Effective          |    Log(T/Kelvin)      |
|                       |    temperature        |                       |
+-----------------------+-----------------------+-----------------------+
|    feh                |    Metallicity        |    [Fe/H]             |
+-----------------------+-----------------------+-----------------------+
|    rad                |    Galactic radial    |    kpc                |
|                       |    distance           |                       |
+-----------------------+-----------------------+-----------------------+
|    isMultiple         |    True if the system |    N/A                |
|                       |    has companions,    |                       |
|                       |    False if the       |                       |
|                       |    system does not    |                       |
+-----------------------+-----------------------+-----------------------+
|    N_companions       |    Number of          |    N/A                |
|                       |    companions         |                       |
+-----------------------+-----------------------+-----------------------+
|    rem_id             |    Integer indicating |    N/A                |
|                       |    the remnant object |                       |
|                       |    type (see          |                       |
|                       |    Additional         |                       |
|                       |    Descriptions       |                       |
|                       |    below)             |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id             |    Object ID-- unique |    N/A                |
|                       |    integer to         |                       |
|                       |    identify           |                       |
|                       |    star/compact       |                       |
|                       |    object             |                       |
+-----------------------+-----------------------+-----------------------+
|    ubv_J, H, K, U, I, |    UBV photometric    |    mag                |
|    B, V, R            |    system, J, H, K,   |                       |
|                       |    U, I, B, V, R      |                       |
|                       |    system absolute    |                       |
|                       |    magnitude          |                       |
+-----------------------+-----------------------+-----------------------+
|    ztf_g, r, i        |    ztf photometric    |    mag                |
|    (optional)         |    system g, r, i     |                       |
|                       |    absoltue magnitude |                       |
+-----------------------+-----------------------+-----------------------+
|    vr                 |    Galactic radial    |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_b               |    Galactic proper    |    mas/yr             |
|                       |    motion, b          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_lcosb           |    Galactic proper    |    mas/yr             |
|                       |    motion, l          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+

..

         Note that the tag names can be used to access HDF5 files (see “How
         to read HDF5 files” above)
         For stars (which are generated by Galaxia ), the following outputs
         are taken directly from Galaxia and just reformatted into the HDF5
         format; parenthetical names correspond to the tag name from Galaxia, 
         if different: zams_mass (smass), mass (mact), px, py, pz, vx, vy,
         vz, age, popid, ubv_k, ubv_i, ubv_u, ubv_b, ubv_v, ubv_r, ubv_j,
         ubv_h, exbv (exbv_schlegel), teff, grav, mbol (lum), feh. Note that
         the lum key from Galaxia is referred to as mbol in the Galaxia
         documentation.
    
        For compact objects (which we generated with our population synthesis
        code, SPISEA ), we must assign these values ourselves.
        
        For both stars and compact objects, the following are things we have
        directly calculated or assigned ourselves: rem_id, rad, glat, glon,
        vr, mu_b, mu_lcosb, obj_id. (For reasons relating to managing RAM, we
        calculate rad, glat, and glon although they are an output given
        directly from Galaxia, and we could have just read in the value.
        However, it can be calculated directly from knowledge of px, py, and
        pz.)

4.1.3 Stars and Compact Objects Companions (HDF5)
++++++++++++++++++++++++++++++++++++++++++++++++++
    
       The data output contained in the HDF5 datasets are a combination of
       outputs that come directly from SPISEA , and outputs we ourselves
       have calculated or defined.

       Default name: *root*\ \_companions.h5

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    system_idx         |    System index       |    N/A                |
|                       |    corresponding to   |                       |
|                       |    the obj_idx of the |                       |
|                       |    primary            |                       |
+-----------------------+-----------------------+-----------------------+
|    zams_mass          |    ZAMS mass          |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    Teff               |    Effective          |    K                  |
|                       |    Temperature        |                       |
+-----------------------+-----------------------+-----------------------+
|    L                  |    Luminosity         |    W                  |
+-----------------------+-----------------------+-----------------------+
|    logg               |    Surface gravity    |    cgs                |
+-----------------------+-----------------------+-----------------------+
|    isWR               |    Is star a          |    N/A                |
|                       |    Wolf-Rayet?        |                       |
+-----------------------+-----------------------+-----------------------+
|    mass               |    Current mass       |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    phase              |    Evolution phase    |    N/A                |
|                       |    (equivalent to     |                       |
|                       |    rem_id in primary  |                       |
|                       |    table)             |                       |
+-----------------------+-----------------------+-----------------------+
|    metallicity        |    Companion          |    [Fe/H]             |
|                       |    metallicity        |                       |
+-----------------------+-----------------------+-----------------------+
|    m_ubv_U, B, V, I,  |    System magnitude   |    mag                |
|    R                  |    in filters from    |                       |
|                       |    SPISEA system      |                       |
+-----------------------+-----------------------+-----------------------+
|    m_ukirt_H, K, J    |    System magnitude   |    mag                |
|                       |    in filters from    |                       |
|                       |    SPISEA system      |                       |
+-----------------------+-----------------------+-----------------------+
|    m_ztf_g, r, i      |    System magnitude   |    mag                |
|                       |    in filters from    |                       |
|                       |    SPISEA system      |                       |
+-----------------------+-----------------------+-----------------------+
|    log_a              |    Log of the system  |    log(AU)            |
|                       |    semimajor axis     |                       |
+-----------------------+-----------------------+-----------------------+
|    e                  |    Eccentricity       |    N/A                |
+-----------------------+-----------------------+-----------------------+
|    i                  |    Inclination        |    deg                |
+-----------------------+-----------------------+-----------------------+
|    Omega              |    Longitude of       |    deg                |
|                       |    ascending node     |                       |
+-----------------------+-----------------------+-----------------------+
|    omega              |    Argument of        |    deg                |
|                       |    periapsis          |                       |
+-----------------------+-----------------------+-----------------------+
|                       |    Difference between |    ΔM⊙                |
|  zams_mass_match_diff |    mass of SPISEA     |                       |
|                       |    primary and        |                       |
|                       |    matched Galaxia    |                       |
|                       |    primary            |                       |
+-----------------------+-----------------------+-----------------------+
|    zams_mass_prim     |    ZAMS mass of       |    M⊙                 |
|                       |    original SPISEA    |                       |
|                       |    priamry            |                       |
+-----------------------+-----------------------+-----------------------+
|    spisea_idx         |    System index in    |    N/A                |
|                       |    original SPISEA    |                       |
|                       |    systems table      |                       |
+-----------------------+-----------------------+-----------------------+

..

4.2 calc_events
----------------

4.2.1 Event candidates table (Astropy FITS table)
+++++++++++++++++++++++++++++++++++++++++++++++++
    
       The event candidates table is very similar to the HDF5 file created
       in perform_pop_syn. (In fact, the top part is completely duplicated;
       it's here for completeness.)
    
       However, the main difference is that there is a LOT less of the
       output, so instead of writing it in arrays in an HDF5 file, we use an
       Astropy table.
    
       Each row in this table is associated with a microlensing event, each
       of which has a lens-source pair
    
       Default name: *root*\ \_events.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    zams_mass (_L,     |    ZAMS mass          |    M⊙                 |
|    \_S)               |                       |                       |
+-----------------------+-----------------------+-----------------------+
|    mass (_L, \_S)     |    Current mass       |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    systemMass (_L,    |    Sum of mass of     |    M⊙                 |
|    \_S)               |    primary and        |                       |
|                       |    companions (if     |                       |
|                       |    existent)          |                       |
+-----------------------+-----------------------+-----------------------+
|    px (_L, \_S)       |    Heliocentric x     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    py (_L, \_S)       |    Heliocentric y     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    pz (_L, \_S)       |    Heliocentric z     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    vx (_L, \_S)       |    Heliocentric x     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vy (_L, \_S)       |    Heliocentric y     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vz (_L, \_S)       |    Heliocentric z     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    age (_L, \_S)      |    Age                |    log(age/yr)        |
+-----------------------+-----------------------+-----------------------+
|    popid (_L, \_S)    |    Population ID -    |    N/A                |
|                       |    integer indicating |                       |
|                       |    the population     |                       |
|                       |    type ranging from  |                       |
|                       |    0 to 9             |                       |
+-----------------------+-----------------------+-----------------------+
|    exbv (_L, \_S)     |    Extinction E(B-V)  |    mag                |
|                       |    at the location of |                       |
|                       |    star given by 3-D  |                       |
|                       |    Schlegel           |                       |
|                       |    extinction maps    |                       |
+-----------------------+-----------------------+-----------------------+
|    glat (_L, \_S)     |    Galactic latitude  |    deg                |
+-----------------------+-----------------------+-----------------------+
|    glon (_L, \_S)     |    Galactic longitude |    deg                |
+-----------------------+-----------------------+-----------------------+
|    mbol (_L, \_S)     |    Bolometric         |    log(L/L⊙)          |
|                       |    magnitude          |                       |
+-----------------------+-----------------------+-----------------------+
|    grav (_L, \_S)     |    Surface gravity    |    log(gravity)       |
+-----------------------+-----------------------+-----------------------+
|    teff (_L, \_S)     |    Effective          |    Log(T/Kelvin)      |
|                       |    temperature        |                       |
+-----------------------+-----------------------+-----------------------+
|    feh (_L, \_S)      |    Metallicity        |    [Fe/H]             |
+-----------------------+-----------------------+-----------------------+
|    rad (_L, \_S)      |    Galactic radial    |    kpc                |
|                       |    distance           |                       |
+-----------------------+-----------------------+-----------------------+
|    isMultiple (_L,    |    True if the system |    N/A                |
|    \_S)               |    has companions,    |                       |
|                       |    False if the       |                       |
|                       |    system does not    |                       |
+-----------------------+-----------------------+-----------------------+
|    N_companions (_L,  |    Number of          |    N/A                |
|    \_S)               |    companions         |                       |
+-----------------------+-----------------------+-----------------------+
|    rem_id (_L, \_S)   |    Integer indicating |    N/A                |
|                       |    the remnant object |                       |
|                       |    type (more details |                       |
|                       |    in tag             |                       |
|                       |    description)       |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id (_L, \_S)   |    Object ID-- unique |    N/A                |
|                       |    integer to         |                       |
|                       |    identify           |                       |
|                       |    star/compact       |                       |
|                       |    object             |                       |
+-----------------------+-----------------------+-----------------------+
|    ubv_J, H, K, U, I, |    UBV photometric    |    mag                |
|    B, V, R (_L, \_S)  |    system, J, H, K,   |                       |
|                       |    U, I, B, V, R      |                       |
|                       |    absolute magnitude |                       |
+-----------------------+-----------------------+-----------------------+
| ztf_g, r, i (_L,      |    ztf photometric    |    mag                |
| \_S)(optional)        |    system g, r, i     |                       |
|                       |    absoltue magnitude |                       |
+-----------------------+-----------------------+-----------------------+
|    vr (_L, \_S)       |    Galactic radial    |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_b (_L, \_S)     |    Galactic proper    |    mas/yr             |
|                       |    motion, b          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_lcosb (_L, \_S) |    Galactic proper    |    mas/yr             |
|                       |    motion, l          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+
|    theta_E            |    (Angular) Einstein |    mas                |
|                       |    radius             |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_rel             |    Relative           |    mas/yr             |
|                       |    source-lens proper |                       |
|                       |    motion             |                       |
+-----------------------+-----------------------+-----------------------+
|    u0                 |    (Unitless) minimum |    | dimensionless    |
|                       |    source-lens        |    | (normalized to   |
|                       |    separation,        |      θE)              |
|                       |    *during* the       |                       |
|                       |    survey             |                       |
+-----------------------+-----------------------+-----------------------+
|    t0                 | Time at which minimum |    days               |
|                       | source-lens           |                       |
|                       | separation occurs     |                       |
+-----------------------+-----------------------+-----------------------+

..

        Tag names ARE used for the Astropy table. You will see a lot of the
        tag names have a parenthetical after (_L, \_S). That is to indicate
        there is one tag for the lens (L) and one for the source (S), since
        for a given event, you need to have both a lens and a source, and
        each of these things has a mass, a velocity, a position, etc. For
        example, zams_mass_L is the ZAMS mass of the lens, and age_S is the
        log(age/yr) of the source.

4.2.2 Blends table (Astropy FITS table)
++++++++++++++++++++++++++++++++++++++++
    
       For each candidate microlensing event, associated with it are blended
       stars, which we call neighbors. Given the blend radius chosen when
       running calc_events, the blend table saves all neighbor stars that
       fall within that distance from the lenses in the candidate events
       table. The blends table is again almost identical to the HDF5 output,
       but is has three additional items. For each neighbor star, it lists
       the object ID of the lens and source it is associated with, and the
       distance between itself and the lens. Note that there can be multiple
       neighbor stars associated with a single lens and source (microlensing
       event).
    
       Default name: *root*\ \_blends.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    zams_mass_N        |    ZAMS mass          |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    mass_N             |    Current mass       |    M⊙                 |
+-----------------------+-----------------------+-----------------------+
|    systemMass_N       |    Sum of mass of     |    M⊙                 |
|                       |    primary and        |                       |
|                       |    companions (if     |                       |
|                       |    existent)          |                       |
+-----------------------+-----------------------+-----------------------+
|    px_N               |    Heliocentric x     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    py_N               |    Heliocentric y     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    pz_N               |    Heliocentric z     |    kpc                |
|                       |    position           |                       |
+-----------------------+-----------------------+-----------------------+
|    vx_N               |    Heliocentric x     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vy_N               |    Heliocentric y     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    vz_N               |    Heliocentric z     |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    age_N              |    Age                |    log(age/yr)        |
+-----------------------+-----------------------+-----------------------+
|    popid_N            |    Population ID -    |    N/A                |
|                       |    integer indicating |                       |
|                       |    the population     |                       |
|                       |    type ranging from  |                       |
|                       |    0 to 9             |                       |
+-----------------------+-----------------------+-----------------------+
|    exbv_N             |    Extinction E(B-V)  |    mag                |
|                       |    at the location of |                       |
|                       |    star given by 3-D  |                       |
|                       |    Schlegel           |                       |
|                       |    extinction maps    |                       |
+-----------------------+-----------------------+-----------------------+
|    glat_N             |    Galactic latitude  |    deg                |
+-----------------------+-----------------------+-----------------------+
|    glon_N             |    Galactic longitude |    deg                |
+-----------------------+-----------------------+-----------------------+
|    mbol_N             |    Bolometric         |    log(L/L⊙)          |
|                       |    magnitude          |                       |
+-----------------------+-----------------------+-----------------------+
|    grav_N             |    Surface gravity    |    log(gravity)       |
+-----------------------+-----------------------+-----------------------+
|    teff_N             |    Effective          |    Log(T/Kelvin)      |
|                       |    temperature        |                       |
+-----------------------+-----------------------+-----------------------+
|    feh_N              |    Metallicity        |    [Fe/H]             |
+-----------------------+-----------------------+-----------------------+
|    rad_N              |    Galactic radial    |    kpc                |
|                       |    distance           |                       |
+-----------------------+-----------------------+-----------------------+
|    isMultiple_N       |    True if the system |    N/A                |
|                       |    has companions,    |                       |
|                       |    False if the       |                       |
|                       |    system does not    |                       |
+-----------------------+-----------------------+-----------------------+
|    N_companions_N     |    Number of          |    N/A                |
|                       |    companions         |                       |
+-----------------------+-----------------------+-----------------------+
|    rem_id_N           |    Integer indicating |    N/A                |
|                       |    the remnant object |                       |
|                       |    type (more details |                       |
|                       |    in tag             |                       |
|                       |    description)       |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_N           |    Object ID-- unique |    N/A                |
|                       |    integer to         |                       |
|                       |    identify           |                       |
|                       |    star/compact       |                       |
|                       |    object             |                       |
+-----------------------+-----------------------+-----------------------+
|    ubv_J, H, K, U, I, |    UBV photometric    |    mag                |
|    B, V, R_N          |    system, J, H, K,   |                       |
|                       |    U, I, B, V, R      |                       |
|                       |    absolute magnitude |                       |
+-----------------------+-----------------------+-----------------------+
|    ztf_g, r, i_N      |    ztf photometric    |    mag                |
|    (optional)         |    system g, r, i     |                       |
|                       |    absoltue magnitude |                       |
+-----------------------+-----------------------+-----------------------+
|    vr_N               |    Galactic radial    |    km/s               |
|                       |    velocity           |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_b_N             |    Galactic proper    |    mas/yr             |
|                       |    motion, b          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+
|    mu_lcosb_N         |    Galactic proper    |    mas/yr             |
|                       |    motion, l          |                       |
|                       |    component          |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_L           |    Object ID of the   |    N/A                |
|                       |    lens               |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_S           |    Object ID of the   |    N/A                |
|                       |    source             |                       |
+-----------------------+-----------------------+-----------------------+
|    sep_LN             |    Separation between |    arcsec             |
|                       |    lens and neighbor  |                       |
+-----------------------+-----------------------+-----------------------+

..

        Note that there is no additional companions table associated with
        calc_events. In order to cross reference between the events and
        
        companions, refine_events must be run first

4.3 refine_events
-----------------

4.3.1 Events table (Astropy FITS table)
++++++++++++++++++++++++++++++++++++++++
    
       The output here is very similar to the candidate events table. In
       fact, part of it is completely duplicated. All tags listed in the
       event candidates table are also part of the events table. However,
       the following columns are also appended. NOTE: the entries for u0 and
       t0 are *overwritten*; the values for u0 and t0 returned from
       calc_events is different from that returned in refine_events. Each
       refine_events file requires you to choose a filter and extinction
       law; in this table we suppose filter *x* is chosen.
    
       Default name: *root*\ \_refine_events\_\ *filter_reddeninglaw*.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag Name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    u0                 |    (Unitless) minimum |    dimensionless      |
|                       |    source-lens        |                       |
|                       |    separation,        |                       |
|                       |    *during* the       |                       |
|                       |    survey             |                       |
+-----------------------+-----------------------+-----------------------+
|    t0                 |    Time at which      |    days               |
|                       |    minimum            |                       |
|                       |    source-lens        |                       |
|                       |    separation occurs  |                       |
+-----------------------+-----------------------+-----------------------+
|    delta_m\_\ *x*     |    Bump amplitude     |    mag                |
|                       |    (difference in     |                       |
|                       |    baseline and       |                       |
|                       |    maximum            |                       |
|                       |    magnification      |                       |
|                       |    magnitude) in      |                       |
|                       |    *x*-band           |                       |
+-----------------------+-----------------------+-----------------------+
|    pi_rel             |    Relative parallax  |    mas                |
+-----------------------+-----------------------+-----------------------+
|    pi_E               |    Microlensing       |    dimensionless      |
|                       |    parallax           |                       |
+-----------------------+-----------------------+-----------------------+
|    t_E                |    Einstein crossing  |    days               |
|                       |    time               |                       |
+-----------------------+-----------------------+-----------------------+
|    ubv\_\ *x*\ \_app  |    UBV photometric    |    mag                |
|    (_L, \_S)          |    system, *x*-band   |                       |
|                       |    apparent           |                       |
|                       |    magnitude, with    |                       |
|                       |    extinction         |                       |
+-----------------------+-----------------------+-----------------------+
|    ubv\_\ *x*\ \_LSN  |    Blended magnitude  |    mag                |
|                       |    in *x*-band        |                       |
|                       |    (Apparent          |                       |
|                       |    magnitude of       |                       |
|                       |    source + lens +    |                       |
|                       |    neighbors →        |                       |
|                       |    “baseline mag”)    |                       |
+-----------------------+-----------------------+-----------------------+
|    f_blend\_\ *x*     |    Source flux        |    dimensionless      |
|                       |    fraction (unlensed |                       |
|                       |    source flux        |                       |
|                       |    divided by         |                       |
|                       |    baseline) in       |                       |
|                       |    *x*-band           |                       |
+-----------------------+-----------------------+-----------------------+
|                       |    Galactic longitude |    deg                |
| cent_glon\_\ *x*\ \_N |    l of neighbor      |                       |
|                       |    stars' centroid    |                       |
+-----------------------+-----------------------+-----------------------+
|                       |    Galactic latitude  |    deg                |
| cent_glat\_\ *x*\ \_N |    l of neighbor      |                       |
|                       |    stars' centroid    |                       |
+-----------------------+-----------------------+-----------------------+
|                       |    Apparent magnitude |    mag                |
|   ubv\_\ *x*\ \_app_N |    of neighbor stars, |                       |
|                       |    *x*-band apparent  |                       |
|                       |    magnitude          |                       |
+-----------------------+-----------------------+-----------------------+
|    pps_seed           |    Seed used in       |    N/A                |
|                       |    perform_pop_syn    |                       |
+-----------------------+-----------------------+-----------------------+
|    gal_seed           |    Seed used in       |    N/A                |
|                       |    run_galaxia        |                       |
+-----------------------+-----------------------+-----------------------+

..

4.3.2 Companions table (Astropy FITS table)
+++++++++++++++++++++++++++++++++++++++++++
    
       This table is very similar to the companion HDF5 file created in
       perform_pop_syn. In fact, part of it is completely duplicated. There
       is some additional information to index between this table and the
       events table and additional binary properties below. There is also no
       mass_match_diff column. Each row in this table is a companion
       associated with an event. So, if a system is lensed twice, its
       companions will be duplicated in this table.
    
       Default name:
       *root*\ \_refine_events\_\ *filter_reddeninglaw\_*\ companions.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    prim_type          |    Type of primary    |    N/A                |
|                       |    associated with    |                       |
|                       |    companion: ‘S' if  |                       |
|                       |    source or 'L’ if   |                       |
|                       |    lens               |                       |
+-----------------------+-----------------------+-----------------------+
|    q                  |    Companion          |    dimensionless      |
|                       |    mass/primary mass  |                       |
+-----------------------+-----------------------+-----------------------+
|    sep                |    Projected angular  |    mas                |
|                       |    separation between |                       |
|                       |    companion and      |                       |
|                       |    primary            |                       |
+-----------------------+-----------------------+-----------------------+
|    P                  |    Period of          |    years              |
|                       |    companion          |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_L           |    Object ID of the   |    N/A                |
|                       |    lens               |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_S           |    Object ID of the   |    N/A                |
|                       |    source             |                       |
+-----------------------+-----------------------+-----------------------+
|    alpha              |    Angle between      |    deg                |
|                       |    binary axis and    |                       |
|                       |    North              |                       |
+-----------------------+-----------------------+-----------------------+
|    phi_pi_E           |    Angle between      |    deg                |
|                       |    North and relative |                       |
|                       |    proper motion      |                       |
|                       |    between the source |                       |
|                       |    and the lens       |                       |
+-----------------------+-----------------------+-----------------------+
|    phi                |    Angle between the  |    deg                |
|                       |    relative proper    |                       |
|                       |    motion and the     |                       |
|                       |    binary axis        |                       |
+-----------------------+-----------------------+-----------------------+

..

4.4 refine_binary_events
-------------------------

   In this section we simulate lightcurves for all the binary events
   that contain a binary and store the parameters. In the case of triple

   lenses/sources we simulate multiple lightcurves and choose the one
   with the largest amplitude. The following are the examples of systems
   we simulate:

    * Binary lens and binary source:
        * Primary lens + companion lens + primary source + companion source
    * Triple lens and single source:
        * Primary lens + companion lens 1 + source
        * Primary lens + companion lens 2 + source
    * Triple lens and binary source:
        * Primary lens + companion lens 1 + primary source + companion source
        * Primary lens + companion lens 2 + primary source + companion source
    * Triple lens and triple source
        * Primary lens + companion lens 1 + primary source + companion source 1
        * Primary lens + companion lens 1 + primary source + companion source 2
        * Primary lens + companion lens 2 + primary source + companion source 1
        * Primary lens + companion lens 2 + primary source + companion source 2

   The parameters for all these lightcurves are stored in the
   lightcurves.fits table (see 4.4.5) where each bullet point would be
   an entry in that file. We then choose the lightcurve with the largest Δm as the
   microlensing event whose parameters are used in other tables. Whether
   it was used or not is indicated in the lightcurve.fits table.

4.4.1 Events table (Astropy FITS table)
++++++++++++++++++++++++++++++++++++++++
    
        This table is a duplicate version of the events table from refine_events
        with some additional properties below from the simulated lightcurves.
    
       Default name:
       *root*\ \_refine_events\_\ *filter_reddeninglaw*\ \_rb.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    n_peaks            |    Number of peaks in |    N/A                |
|                       |    lightcurve         |                       |
+-----------------------+-----------------------+-----------------------+
|    bin_delta_m        |    Bump amplitude     |    mag                |
|                       |    (difference in     |                       |
|                       |    baseline and       |                       |
|                       |    maximum            |                       |
|                       |    magnification      |                       |
|                       |    magnitude)         |                       |
+-----------------------+-----------------------+-----------------------+
|    tE_sys             |    | Empirical        |    days               |
|                       |      Einstein         |                       |
|                       |      crossing time    |                       |
|                       |      (when the system |                       |
|                       |      magnitude is at  |                       |
|                       |      least 10% the    |                       |
|                       |    | maximum          |                       |
|                       |      magnitude)       |                       |
+-----------------------+-----------------------+-----------------------+
|    tE_primary         |    Empirical Einstein |    days               |
|                       |    crossing time of   |                       |
|                       |    the peak of max    |                       |
|                       |    mag (when the      |                       |
|                       |    system magnitude   |                       |
|                       |    is at least 50%    |                       |
|                       |    the maximum        |                       |
|                       |    magnitude of peak) |                       |
+-----------------------+-----------------------+-----------------------+
|    primary_t          |    Time at which      |    days               |
|                       |    maximum peak       |                       |
|                       |    occurs             |                       |
+-----------------------+-----------------------+-----------------------+
|    avg_t              |    Average time the   |    days               |
|                       |    peaks occur        |                       |
+-----------------------+-----------------------+-----------------------+
|    std_t              |    Standard deviation |    days               |
|                       |    of times peaks     |                       |
|                       |    occur              |                       |
+-----------------------+-----------------------+-----------------------+
|    asymmetry          |    Asymmetry as       |    dimensionless      |
|                       |    defined by         |                       |
|                       |    Chebyshev          |                       |
|                       |    Polynomials (see   |                       |
|                       |    Additional         |                       |
|                       |    Descriptions       |                       |
|                       |    below)             |                       |
+-----------------------+-----------------------+-----------------------+
|    companion_idx_list |    List of companion  |    N/A                |
|                       |    indices associated |                       |
|                       |    with events        |                       |
|                       |    (corresponds with  |                       |
|                       |    companion_idx in   |                       |
|                       |    the companions     |                       |
|                       |    table)             |                       |
+-----------------------+-----------------------+-----------------------+

..

4.4.2 Companions table (Astropy FITS table)
+++++++++++++++++++++++++++++++++++++++++++++
         This table is a duplicate version of the companions table from
         refine_events with an additional id Default name:
         *root*\ \_refine_events\_\ *filter_reddeninglaw\_*\ companions_rb.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    companion_idx      |    Companion index    |    N/A                |
|                       |    which corresponds  |                       |
|                       |    to position in the |                       |
|                       |    array              |                       |
+-----------------------+-----------------------+-----------------------+

..

4.4.3 Multi-peak table (Astropy FITS table)
++++++++++++++++++++++++++++++++++++++++++++
         This table describes properties of each of the peaks in a binary
         lens microlensing event lightcurve that passed a significance
         threshold of 10-5 (by default).

       Default name:
       *root*\ \_refine_events\_\ *filter_reddeninglaw\_*\ companions_rb_mp.fits

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    comp_id            |    Companion ID - the |    N/A                |
|                       |    position of        |                       |
|                       |    associated         |                       |
|                       |    companion in       |                       |
|                       |    companion table    |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_L           |    Object ID of the   |    N/A                |
|                       |    lens               |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_S           |    Object ID of the   |    N/A                |
|                       |    source             |                       |
+-----------------------+-----------------------+-----------------------+
|    n_peaks            |    Number of peaks in |    N/A                |
|                       |    lightcurve         |                       |
+-----------------------+-----------------------+-----------------------+
|    t                  |    Time at which peak |    days               |
|                       |    occurs             |                       |
+-----------------------+-----------------------+-----------------------+
|    tE                 |    Empirical Einstein |    days               |
|                       |    crossing time of   |                       |
|                       |    the peak (when the |                       |
|                       |    system magnitude   |                       |
|                       |    is at least 50%    |                       |
|                       |    the maximum        |                       |
|                       |    magnitude of peak) |                       |
+-----------------------+-----------------------+-----------------------+
|    delta_m            | Bump amplitude        |    mag                |
|                       | (difference in        |                       |
|                       | baseline and maximum  |                       |
|                       | magnification         |                       |
|                       | magnitude of peak)    |                       |
+-----------------------+-----------------------+-----------------------+
|    ratio              |    Magnitude ratio    |    dimensionless      |
|                       |    between peak of    |                       |
|                       |    the maximum mag    |                       |
|                       |    and this peak      |                       |
+-----------------------+-----------------------+-----------------------+

..

4.4.4 Lightcurve table (Astropy FITS table)
++++++++++++++++++++++++++++++++++++++++++++
    
        This table has the parameters of all the lightcurves simulated. In any
        cases involving a binary lens or source but no triples, there will be
        only one lightcurve for that event. However, in triple lens or triple
        source cases there will be 2 lightcurves in the TSBL, TSPL, BSTL, and
        PSTL cases and 4 lightcurves in the TSTL case.

+-----------------------+-----------------------+-----------------------+
|    **Tag name**       |    **Brief            |    **Units**          |
|                       |    Description**      |                       |
+=======================+=======================+=======================+
|    obj_id_L           |    Object ID of the   |    N/A                |
|                       |    lens               |                       |
+-----------------------+-----------------------+-----------------------+
|    obj_id_S           |    Object ID of the   |    N/A                |
|                       |    source             |                       |
+-----------------------+-----------------------+-----------------------+
|    companion_L        |    Companion ID of    |    N/A                |
|                       |    the lens (blank if |                       |
|                       |    BSPL)              |                       |
+-----------------------+-----------------------+-----------------------+
|    companion_S        | Companion ID of the   |    N/A                |
|                       | source (blank if      |                       |
|                       | PSBL)                 |                       |
+-----------------------+-----------------------+-----------------------+
|    class              |    Type of            |    N/A                |
|                       |    microlensing event |                       |
|                       |    simulated (PSBL,   |                       |
|                       |    BSPL, or BSBL)     |                       |
+-----------------------+-----------------------+-----------------------+
|    n_peaks            |    Number of peaks in |    N/A                |
|                       |    lightcurve         |                       |
+-----------------------+-----------------------+-----------------------+
|    bin_delta_m        |    Bump amplitude     |    mag                |
|                       |    (difference in     |                       |
|                       |    baseline and       |                       |
|                       |    maximum            |                       |
|                       |    magnification      |                       |
|                       |    magnitude)         |                       |
+-----------------------+-----------------------+-----------------------+
|    tE_sys             |    | Empirical        |    days               |
|                       |      Einstein         |                       |
|                       |      crossing time    |                       |
|                       |      (when the system |                       |
|                       |      magnitude is at  |                       |
|                       |      least 10% the    |                       |
|                       |    | maximum          |                       |
|                       |      magnitude)       |                       |
+-----------------------+-----------------------+-----------------------+
|    tE_primary         |    Empirical Einstein |    days               |
|                       |    crossing time of   |                       |
|                       |    the peak of max    |                       |
|                       |    mag (when the      |                       |
|                       |    system magnitude   |                       |
|                       |    is at least 50%    |                       |
|                       |    the maximum        |                       |
|                       |    magnitude of peak) |                       |
+-----------------------+-----------------------+-----------------------+
|    primary_t          |    Time at which      |    days               |
|                       |    maximum peak       |                       |
|                       |    occurs             |                       |
+-----------------------+-----------------------+-----------------------+
|    avg_t              |    Average time the   |    days               |
|                       |    peaks occur        |                       |
+-----------------------+-----------------------+-----------------------+
|    std_t              |    Standard deviation |    days               |
|                       |    of times peaks     |                       |
|                       |    occur              |                       |
+-----------------------+-----------------------+-----------------------+
|    asymmetry          |    Asymmetry as       |    dimensionless      |
|                       |    defined by         |                       |
|                       |    Chebyshev          |                       |
|                       |    Polynomials (see   |                       |
|                       |    Additional         |                       |
|                       |    Descriptions       |                       |
|                       |    below)             |                       |
+-----------------------+-----------------------+-----------------------+
|    used_lightcurve    |    If this lightcurve |    N/A (bool)         |
|                       |    had the largest    |                       |
|                       |    bin_delta_m of the |                       |
|                       |    set, it will be    |                       |
|                       |    used in the other  |                       |
|                       |    tables             |                       |
|                       |    corresponding to   |                       |
|                       |    this event (True). |                       |
|                       |    If not, False.     |                       |
+-----------------------+-----------------------+-----------------------+

4.5 Additional descriptions of tags
------------------------------------

    **rem_id/phase:** These label the different types of remnant objects (star, black hole, neutron star, or white dwarf). 
    They are identified
    as following:
    
        * 0: Star
        * 101: White Dwarf
        * 102: Neutron Star
        * 103: Black Hole
    
    **pop_id:** Describes which population is generated.
    
        * 0: Thin disk, ≤0.15 Gyr
        * 1: Thin disk, 0.15-1 Gyr
        * 2: Thin disk, 1-2 Gyr
        * 3: Thin disk, 2-3 Gyr
        * 4: Thin disk, 3-5 Gyr
        * 5: Thin disk, 5-7 Gyr
        * 6: Thin disk, 7-10 Gyr
        * 7: Thick disk, 11 Gyr (single-age)
        * 8: Stellar halo, 14 Gyr (single-age)
        * 9: Bulge, 10 Gyr (single-age)
        * In Galaxia there is an option for a 10th population type; the Bullock and Johnston stellar halos. We have chosen not use it, and the code is not written to include it.
    
    **px, py, pz; vx, vy, vz:** These are given in heliocentric coordinates
    (i.e. Cartesian coordinates with the sun at the origin.) See subsection
    on coordinate systems for more information.
    
    **rad, glat, glon; vr, mu_b, mu_lcosb:** These are given in galactic
    coordinates (i.e. spherical coordinates with the sun at the origin.) See
    subsection on coordinate systems for more information.
    
    **ubv_U, B, V, R, I, J, H, K; exbv:** Photometry information is given in
    absolute magnitude. For NSs and BHs, all these values are \\texttt{nan}
    to indicate they are dark. Note that for dark primaries with luminous
    companions, these values will be the system luminosity.
    
    **t0:** Note that you can have a negative day (this just means time
    before the “zero" time, which is defined as the state of the system that
    is generated by Galaxia and the population synthesis. Since we are
    assuming everything moves in straight lines, we can propagate either
    forward or backwards.) This can also be the case for **primary_t** and
    **t** in the output of refine_binary_events.
    
    **asymmetry:** Asymmetry as defined by (`Night et al. 2008 <https://iopscience.iop.org/article/10.1086/590320>`_) where the light curve is fit
    by a Chebyshev polynomial to the 50th degree (`Khakpash et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021AJ....161..132K/abstract>`_) and :math:`k = \left(\sum T^2_{\rm odd} / \sum T^2_{\rm even} \right)`.
    Here the T's are the coefficients of the Chebyshev polynomial. 
    A light curve is symmetric when k = 0.
    
====================
5 Coordinate Systems
====================
There are two different coordinate systems used, Heliocentric and Galactic. 
Heliocentric coordinates are Cartesian coordinates with the sun at the origin. The positive :math:`x` axis is pointing toward the Galactic Center, and the positive :math:`z` axis is pointing toward the Galactic North Pole.
Galactic coordinates are spherical coordinates with the sun at the origin. 
Longitude :math:`l` is measuring the angular distance of an object eastward along the galactic equator from the galactic center, and latitude :math:`b` is measuring the angle of an object north or south of the galactic equator (or midplane) as viewed from Earth; positive to the north, negative to the south. 
Radius :math:`r` is the distance from the sun to the object. 

The conversion between Heliocentric and Galactic is just the same as converting between rectangular to spherical coordinates, where :math:`\phi = l` and :math:`\theta = -b + 90^{\circ}`.
Going from Galactic to Heliocentric (units are degrees):

    :math:`x = r\sin(-b + 90^{\circ}) \cos l =  r \cos b \cos l`

    :math:`y = r\sin(-b + 90^{\circ}) \sin l = r \cos b \sin l`

    :math:`z = r\cos(-b + 90^{\circ}) = r\sin b`

Going from Heliocentric to Galactic (units are degrees):

    :math:`r = \sqrt{x^2 + y^2 + z^2}`

    :math:`b = -\cos^{-1}(z/r) + 90^{\circ}`

    :math:`l = \tan^{-1}(y/x)`

Note: be careful with the branch of arctangent. Practically, use ``numpy.arctan2`` if using Python.

.. image:: popsycle_docs_images/media/coords.png
   :align: center

Diagram of Heliocentric and Galactic coordinate systems. The red dot is the sun.

