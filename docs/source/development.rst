.. _development:

Development Process for PopSyCLE
================================

We are excited to welcome anyone who wants to develop new features for PopSyCLE. In order to keep
this process organized, we ask that you follow this process for creating and submitting new features
to the repository.

1 Getting Repository Access
----------------------------
Your first step will be to get added to the repository as a collaborator. Email jlu.astro@berkeley.edu
with the subject line Requesting PopSyCLE Repository Access and include the following information
in your request:

    #. Name
    #. GitHub Username
    #. Academic Affiliation
    #. Summary of the contribution you wish to make

2 Creating a New Feature
------------------------
Due to the small size of the current PopSyCLE community we have opted to develop new features using
branches on this repository instead of using GitHub forks. Begin by 
`cloning the current repository <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_
onto your local machine and checking out the main branch. From main, create a new branch with
a succinct name summarizing the feature you plan to write. Push this new branch back to GitHub
so that all collaborators can see that a new feature is in development.

3 Submitting a New Feature for Review
--------------------------------------
Once you have completed the development of a new feature, you will want to submit it for review to
be added into the master branch. We use GitHub pull requests for such review. `Create a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_
on GitHub asking to merge your feature branch back into the master branch. Your feature will then
be examined by one of our approved reviewers and you may be asked to comment on a `pull request
review <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#about-reviewing-pull-requests>`_. Once the reviewers have determined that your feature is complete without issues, your pull
request will be approved and your branch will be merged into master. Approved reviewers are:

    #. Jessica Lu : `jluastro <https://github.com/jluastro>`_
    #. Casey Lam : `caseylam <https://github.com/caseylam>`_
    #. Natasha Abrams : `nsabrams <https://github.com/nsabrams>`_

You are encouraged to `add these GitHub accounts as reviewers onto your pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/requesting-a-pull-request-review>`_.

4 Reviewing a Pull Request
--------------------------
Pull requests are evaluated to make sure that they successfully implement the intended feature without
breaking any existing code. You can run the tests in the tests directory to verify this.
Reviewers should be mindful to think through all possible consequences
of the proposed code changes, even to files not included in the pull request. If an issue is discovered,
reviewers should `start a review <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#starting-a-review>`_ including commenting on specific lines of code that have issues. Work
with the author of the changes to get the code to a stable and useful state. Then `approve any reviews <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/approving-a-pull-request-with-required-reviews>`_
and `merge the pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request>`_. Pull requests should be merged with the create a merge commit option
and the feature branch should be tagged and archived.
