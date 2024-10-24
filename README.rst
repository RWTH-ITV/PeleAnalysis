
PeleAnalysis
============

This repository contains a collection of standalone routines for processing plotfiles created with the AMReX software framework for block-structured adaptive mesh refinement simulations.  Documentation is under development, but is
available at ::

   https://peleanalysis.readthedocs.io/en/latest/

AMReX is required for these tools, and is available https://github.com/AMReX-Codes/amrex

To clone this repository including the required submodules such as PelePhysics, AMReX, and qslim, use the ``--recursive`` option of the git clone command. The additional options ``--shallow-submodules`` and ``--single-branch`` help decreasing the size of your clone. ::

   git clone --recursive --shallow-submodules --single-branch git@github.com:ITV-RWTH/PeleAnalysis.git

Using this command will automatically like your PeleAnalysis Repository to the repositories of the submodules. Alternatively, you can include your version of PelePysics and AMReX used in your Pele installation by setting the variables ``PELE_PHYSICS_HOME`` and ``AMREX_HOME`` in your ``GNUmakefile`` to point to the local folder where that is placed. After cloning this repository, ``cd Src`` and edit the ``GNUmakefile`` to select which tool to build.  If AMReX is configured properly, a stand-alone executable will be built locally, based on the selected options, including spatial dimension (2 or 3), compiler choices, whether to build with MPI and/or OpenMP enabled, and whether to build a debugging or optimized version.  Note that some of the tools require building a companion f90 source file - you must manually set the flag in the ``GNUmakefile`` accordingly.  More extensive documentation is available (see building instructions below).

Contributions
-------------

To add a new feature to PeleAnalysis, the procedure is:

1. Create a branch for the new feature (locally) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the ``development`` branch into your ``AmazingNewFeature`` branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleAnalysis repository (if you have write access, otherwise fork the repo and
push the new branch to your fork)::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4.  Submit a merge request through the github project page - be sure you are requesting to merge your branch to the ``development`` branch of the ``ITV-RWTH/PeleAnalysis`` and not to the main repository ``AMReX-Combustion/PeleAnalysis``.




Documentation
-------------
Documentation for the analysis routines exists in the Docs directory. To build the documentation::

    cd Docs
    make html


Acknowledgment
--------------
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
