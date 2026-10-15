.. _development:

Development Process for PopSyCLE
================================

We are excited to welcome anyone who wants to develop new features for PopSyCLE. In order to keep
this process organized, we ask that you follow this process for creating and submitting new features
to the repository.

1. Getting Repository Access
----------------------------
Your first step will be to get added to the repository as a collaborator. Email jlu.astro@berkeley.edu
with the subject line Requesting PopSyCLE Repository Access and include the following information
in your request:

    #. Name
    #. GitHub Username
    #. Academic Affiliation
    #. Summary of the contribution you wish to make

2. Creating a New Feature
-------------------------
Due to the small size of the current PopSyCLE community we have opted to develop new features using
branches on this repository instead of using GitHub forks. Begin by 
`cloning the current repository <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_
onto your local machine and checking out the main branch. From dev, create a new branch with
a succinct name summarizing the feature you plan to write. Push this new branch back to GitHub
so that all collaborators can see that a new feature is in development. 
We also ask that you indicate that you are working on this feature in `GitHub issues <https://github.com/jluastro/PopSyCLE/issues>`_.

3. Submitting a New Feature for Review
--------------------------------------
Once you have completed the development of a new feature, you will want to submit it for review to
be added into the dev branch, and then the main branch. We use GitHub pull requests for such review. `Create a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_
on GitHub asking to merge your feature branch back into the master branch. Your feature will then
be examined by one of our approved reviewers and you may be asked to comment on a `pull request
review <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#about-reviewing-pull-requests>`_. Once the reviewers have determined that your feature is complete without issues, your pull
request will be approved and your branch will be merged into master. Approved reviewers are:

    #. Jessica Lu : `jluastro <https://github.com/jluastro>`_
    #. Casey Lam : `caseylam <https://github.com/caseylam>`_
    #. Natasha Abrams : `nsabrams <https://github.com/nsabrams>`_

You are encouraged to `add these GitHub accounts as reviewers onto your pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/requesting-a-pull-request-review>`_.

4. Reviewing a Pull Request
---------------------------
Pull requests are evaluated to make sure that they successfully implement the intended feature without
breaking any existing code. You can run the tests in the tests directory to verify this.
Reviewers should be mindful to think through all possible consequences
of the proposed code changes, even to files not included in the pull request. If an issue is discovered,
reviewers should `start a review <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#starting-a-review>`_ including commenting on specific lines of code that have issues. Work
with the author of the changes to get the code to a stable and useful state. Then `approve any reviews <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/approving-a-pull-request-with-required-reviews>`_
and `merge the pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request>`_. Pull requests should be merged with the create a merge commit option
and the feature branch should be tagged and archived.

5. Building PopsyCLE
--------------------
First, make sure your feature branch is successfully merged into the ``dev``
branch through an approved pull request. If this is your first time building PopSyCLE,
please see the `First Time Conda Setup`_ section as well.

**Increment the version number.**

Edit the version number in the following two files.
Note you MUST increase the number every time you push to pip. This example uses
version number of 1.0.1 -- you will need to change this in the following files:

- popsycle/popsycle/__init__.py
- popsycle/pyproject.toml

Be sure to commit these changes into dev and then merge to main.

**Build the distributions**

Make the wheels and tar.gz files files from the ``main`` branch:

.. code-block:: shell

    cd popsycle
    python -m build

**Prepare the upload to PyPI**

You will need to have a PyPI account (use VIP for authentication
token) and an API token per instructions
`here <https://packaging.python.org/en/latest/tutorials/packaging-projects/>`_.
Note, check the ``.pypirc`` file in your home directory if you have done this before.

**Upload the distribution to PyPI**

Upload the ``.tar.gz`` file to PyPI with the following command:

.. code-block:: shell

   twine upload --repository pypi dist/popsycle-1.0.1.*

The binary should be available to pip install now.

**Update the conda recipe**

In order to update the conda recipe, you will need to increment the version number and hash.
First, pull the SHA256 key from PyPI JSON:

.. code-block:: shell

   curl -s https://pypi.org/pypi/popsycle/1.0.1/json | python3 -c "import sys, json; print(json.load(sys.stdin)['urls'][0]['digests']['sha256'])"

Then, go to the conda recipe:

.. code-block:: shell

   cd ~/pycode/staged-recipes/recipes/popsycle/

Edit the ``meta.yaml`` file to up the version number
(matched to the new pip version number) and hash.
Commit and push the recipe changes up to conda-forge.

.. code-block:: shell

   git commit -a -m "Updating to version 1.0.1"
   git push

**Make a release on GitHub**
On the GitHub PopSyCLE page, do the following:
- Click on Releases on the the left-hand bar.
- Click the Draft a New Release button.
- Make a new tag with something like “v1.0.5”.
- Fill in the description.
- Make the release.

This should automatically be uploaded to Zenodo and assigned a DOI.

First Time Conda Setup
^^^^^^^^^^^^^^^^^^^^^^
You will have to follow the instructions
`here <https://www.pyopensci.org/python-package-guide/tutorials/publish-conda-forge.html>`_
to fork the conda-forge staged-recipes repository, clone it locally,
and make your own branch called popsycle.

Fork and clone https://github.com/conda-forge/staged-recipes
and make a new branch called ``popsycle`` from main:

.. code-block:: shell

   git branch popsycle
   git checkout popsycle

Prepare the popsycle recipe:

.. code-block:: shell

   cd staged-recipes
   mkdir popsycle
   cd popsycle
   echo "" > recipe.yaml

Modify the recipe.yaml (see example `here <https://conda-forge.org/docs/maintainer/adding_pkgs/#id4>`_).
Pull the SHA256 key from PyPI JSON:

.. code-block:: shell

   curl -s https://pypi.org/pypi/popsycle/1.0.1/json | python3 -c "import sys, json; print(json.load(sys.stdin)['urls'][0]['digests']['sha256'])"

Make sure to update the version number.


Build Resources
^^^^^^^^^^^^^^^
- https://packaging.python.org/en/latest/tutorials/packaging-projects/
- https://www.pyopensci.org/python-package-guide/tutorials/publish-conda-forge.html
- https://conda-forge.org/docs/maintainer/adding_pkgs/#step-by-step-instructions


