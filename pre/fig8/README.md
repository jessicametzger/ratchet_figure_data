# Revisiting the Ratchet Principle (companion), fig 8 #

This directory contains code and data necessary for generating Fig. 8.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash scripts ``data/Tofx/loop.sh`` (spatially-varying temperature) and ``data/T_const/loop.sh`` (constant temperature). It generates data that it saves in the respective folders. A copy of the data is already there.

The figure is generated using the script ``fig8.py``, which uses the data in the ``data/T_const/`` and ``data/Tofx/`` folders.

* ``fig8.py`` Code to create figure
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
    * ``pbp-ni-Fofx-disp-prof-J-F-T `` Compiled interacting PBP simulator without spatially-varying temperature
    * ``pbp-ni-Tofx-Fofx-disp-prof-J-F-T`` Compiled non-interacting PBP simulator with spatially-varying temperature
* ``data/``
    * ``Tofx/``
        * ``run.sh`` script that runs ``pbp-ni-Tofx-Fofx-disp-prof-J-F-T`` with provided parameters
        * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-Jprof`` thermal/active stress density data for seed ``XXXXX``
        * ``XXXXX-Fprof`` force density data for seed ``XXXXX``
        * ``XXXXX-Tprof`` average temperature profile for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
    * ``T_const/``
        * ``run.sh`` script that runs ``pbp-ni-Fofx-disp-prof-J-F-T`` with provided parameters
        * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-Jprof`` thermal/active stress density data for seed ``XXXXX``
        * ``XXXXX-Fprof`` force density data for seed ``XXXXX``
        * ``XXXXX-Tprof`` average temperature profile for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
