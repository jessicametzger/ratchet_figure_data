# Revisiting the Ratchet Principle (companion), fig 10 #

This directory contains code and data necessary for generating Fig. 10.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash scripts ``data/asymm_T-step_V/loop.sh`` (asymmetric T(x) and symmetric V(x)), ``data/step_T-symm_V/loop.sh`` (symmetric V(x) and symmetric step T(x) with mutual phase shift), and ``data/symm_T-symm_V/loop.sh`` (symmetric oscillating T(x) and V(x) with mutual phase shift). It generates data that it saves in the respective folders. A copy of the data is already there.

The figure is generated using the script ``fig10.py``, which uses the data in the ``data/asymm_T-step_V/``, ``data/symm_T-symm_V/``, and ``data/step_T-symm_V/`` folders.

* ``fig10.py`` Code to create figure
* ``C_code/``
    * ``util/``
        * ``inhom_fluct.c`` functions to define cubic & linear inhomogeneous temperature & force landscapes
        * ``PFs.c`` functions for calculating pairwise forces
        * ``spatial_hashing.c`` functions for looping through particle force calculations using spatial hashing algorithm
        * ``util.c`` miscellaneous functions
        * ``mt64.h`` For the random-number generator
        * ``mt19937-64.c`` For the random-number generator
    * ``OD-PBPs-2D.c`` Main simulation code for PBPs
    * ``OD-PBPs-2D-functions.c`` Functions for the PBP simulations
    * ``options.h`` Pre-compiler options for the PBP simulations
    * ``pbp-ni-Tofx-Fofx-disp-prof-J-F-T`` Compiled non-interacting PBP simulator with spatially-varying temperature and forces
* ``data/``
    * ``asymm_T-step_V/``
        * ``run.sh`` script that runs ``pbp-ni-Tofx-Fofx-disp-prof-J-F-T`` with provided parameters
        * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-Jprof`` thermal/active stress density data for seed ``XXXXX``
        * ``XXXXX-Fprof`` force density data for seed ``XXXXX``
        * ``XXXXX-Tprof`` average temperature profile for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
    * ``step_T-symm_V/``
        * same as above
    * ``symm_T-symm_V/``
        * same as above