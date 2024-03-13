PopSyCLE
--------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

**Pop**\ulation **Sy**\nthesis for **C**\ompact-object **L**\ensing **E**\vents
is a code to simulate a model of the Milky Way including compact objects and multiple systems and perform a mock microlensing survey.                                               
You can use it to put black hole candidates into context and to understand the effect of Galactic properties on photometric and astrometric microlensing simulation distributions among many other applications.

Full PopSyCLE documentation can be found `here <https://popsycle.readthedocs.io/en/latest/>`_.

Dependencies
------------
`galaxia <http://galaxia.sourceforge.net>`_
PopSyCLE requires a custom version of galaxia in order to support
user selected galaxy models. Please follow the installation instructions
found at our galaxia GitHub repo: https://github.com/jluastro/galaxia.

`SPISEA <https://spisea.readthedocs.io/en/latest/>`_

`BAGLE <https://github.com/MovingUniverseLab/BAGLE_Microlensing>`_


Installation
------------

To install PopSyCLE, clone the GitHub repository and add the repository's
path to your `PYTHONPATH`. For example:

.. code-block:: bash

    git clone git@github.com:jluastro/PopSyCLE.git
    echo "export PYTHONPATH=$PWD/PopSyCLE:$PYTHONPATH" >> ~/.bashrc

Running PopSyCLE
----------------

An example of implementing PopSyCLE can be found
in our example notebook for `single stars <docs/PopSyCLE_example.ipynb>`_ or 
`with multiples <docs/PopSyCLE_example_multiples.ipynb>`_.

Running the PopSyCLE Pipeline
-----------------------------

An example of running the PopSyCLE pipline from the command line be found
`in our example notebook <docs/PopSyCLE_example_run.ipynb>`_.

Running the PopSyCLE Pipeline on a Slurm Scheduler
--------------------------------------------------

An example of running the PopSyCLE pipline on a compute cluster with a
slurm scheduler can be found
`in our example notebook <docs/PopSyCLE_example_slurm.ipynb>`_.

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

`Check out all of PopSyCLE's contributors! <contributors.md>`_.

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

We have outlined the process for developing new features, located at
`<docs/Documentation__PopSyCLE_Development.pdf>`_.
