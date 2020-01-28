PopSyCLE
--------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

**Pop**\ulation **Sy**\nthesis for **C**\ompact-object **L**\ensing **E**\vents


Dependencies
------------
`Galaxia <http://galaxia.sourceforge.net>`_

`PyPopStar <https://pypopstar.readthedocs.io/en/latest/>`_


Installation
------------

PopSyCLE can be installed as a standard python package.

Example
-------

An example of implementing PopSyCLE can be found
`in our example notebook <docs/PopSyCLE_example.ipynb>`_.

Running PopSyCLE on a Slurm Scheduler
-------------------------------------

Slurm is an open source, fault-tolerant, and highly scalable cluster management
and job scheduling system for large and small Linux clusters
(https://slurm.schedmd.com/overview.html). Multiple instances of the PopSyCLE
pipeline can be executed in parallel if running on a system that both has
PopSyCLE and a slurm scheduler installed. The
``popsycle.slurm.generate_slurm_scripts`` function will create and submit
slurm batch scripts that execute the entire PopSyCLE pipeline.

To begin, create a slurm configuration file that contains the linux cluster
parameters necessary to submit a batch scripts. These parameters outline
features of the slurm scheduler and the config file should only need to be
created once. An example can be
found in `<popsycle/data/slurm_config.yaml>`_ and a custom slurm config file
can be generated using ``popsycle.synthetic.generate_slurm_config_file``.
Next, create a PopSyCLE
configuration file that contains the parameters necessary to run the PopSyCLE
pipeline. You may want to create a single PopSyCLE configuration file for
each project in order to keep the execution of PopSyCLE the same over multiple
fields. An example can be found in
`<popsycle/data/popsycle_config.yaml>`_ and a custom popsycle config file
can be generated using ``popsycle.synthetic.generate_popsycle_config_file``.
Lastly, run the
``popsycle.synethtic.generate_slurm_scripts`` function along with parameters
about a single PopSyCLE field that you would like to run. These parameters
could be different for each field and allow for the user to launch multiple
instances of PopSyCLE at once by running this function multiple times with
different input.

Running ``generate_slurm_scripts`` creates a slurm batch script that is
submitted to the slurm scheduler by deafult.

The function ``generate_slurm_scripts`` can be executed by running:

.. code-block:: python

   from popsycle.synthetic import generate_slurm_scripts
   generate_slurm_scripts(slurm_config_filename='slurm_config.yaml',
                          popsycle_config_filename='popsycle_config.yaml',
                          longitude=10.0,
                          latitude=5.0,
                          area=0.33,
                          n_cores_calc_events=4,
                          walltime='01:00:00,
                          seed=None, overwrite=True, submitFlag=True):

The ``generate_slurm_scripts`` docstring explains the function's parameters::

    Generates the slurm script that executes the PopSyCLE pipeline

    Parameters
    ----------
    slurm_config_filename : str
        Name of slurm_config.yaml file containing the slurm parameters
        that will be used the generate the slurm script header.

    popsycle_config_filename : str
        Name of popsycle_config.yaml file containing the PopSyCLE parameters
        that will be passed along to the run_on_slurm.py command in the
        slurm script.

    path_run : str
        Directory containing the parameter file and PopSyCLE output files

    output_root : str
        Base filename of the output files
        Examples:
           '{output_root}.h5'
           '{output_root}.ebf'
           '{output_root}_events.h5'

    longitude : float
        Galactic longitude, ranging from -180 degrees to 180 degrees

    latitude : float
        Galactic latitude, ranging from -90 degrees to 90 degrees

    area : float
        Area of the sky that will be generated, in square degrees

    n_cores_calc_events : int
        Number of cores for executing synthetic.calc_events

    walltime : str
        Amount of walltime that the script will request from slurm.

    Optional Parameters
    -------------------
    seed : int
        If set to non-zero, removes all random sampling and forces identical
        output for Galaxia, PyPopStar and PopSyCLE.
        Default None.

    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    submitFlag : bool
        If set to True, script will be submitted to the slurm scheduler
        after being written to disk. If set to False, it will not be submitted.
        Default is True

The PopSyCLE pipeline can also be run without using slurm scripts using the
executable located at `<popsycle/run.py>`_. Running this executable from the
command line requires the creation of a field configuration file using
```popsycle.synthetic.generate_field_config_file```. More details on
how to run this executable can be found by running
```python {PATH_TO_POPSYCLE}/popsycle/run.py -h```.

License
-------

This project is Copyright (c) Casey Lam and Jessica Lu and licensed under
the terms of the GNU GPL v3+ license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.


Contributing
------------

We love contributions! PopSyCLE is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
popsycle based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
