# Revisiting the Ratchet Principle (companion), fig 3 #

This directory contains code and data necessary for generating Fig. 3.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash script ``data/loop.sh``. It generates data that it saves in the folder ``data/``. A copy of the data is already there.

The data is parsed using the script ``data/get_data.py`` and left in the folders ``data/int/`` and ``data/nonint/``. A copy of the parsed data is already there.

Finally, the figure is generated using the script ``fig3.py``, which uses the parsed data in the ``data/`` folder.

* ``fig3.py`` Code to create figure
* ``C_code/``
    * ``util/``
        * ``inhom_fluct.c`` functions to define cubic inhomogeneous fluctuation landscape
        * ``PFs.c`` functions for calculating pairwise forces
        * ``spatial_hashing.c`` functions for looping through particle force calculations using spatial hashing algorithm
        * ``util.c`` miscellaneous functions
        * ``mt64.h`` For the random-number generator
        * ``mt19937-64.c`` For the random-number generator
    * ``OD-PBPs-2D.c`` Main simulation code for PBPs
    * ``OD-PBPs-2D-functions.c`` Functions for the PBP simulations
    * ``options.h`` Pre-compiler options for the PBP simulations
    * ``pbp-Tofx-Uofx`` Compiled interacting PBP simulator
* ``data/``
    * ``get_data.py`` Gets the results from the C simulations and outputs in the format needed for the ``fig3.py`` script
    * ``run.sh`` script that runs simulation of interacting PBPs with provided discretization and seed. To implement the discretization, a potential landscape is created such that $U'(x) = -\alpha T'(x)$.
    * ``loop.sh`` script that repeatedly calls ``run.sh`` for different discretizations and random number seeds
    * ``0.0/`` discretization $\alpha=0$ data
        * ``disp_avg`` averaged particle displacements for evenly-spaced time intervals
        * ``disp_err`` standard error of the mean average particle displacement over all seeds, for evenly-spaced time intervals
        * ``prof_avg`` average particle density for evenly spaced $x$ from 0 to Lx
        * ``XXXXXX-disp`` raw particle displacement data for seed ``XXXXXX``
        * ``XXXXXX-pos`` raw particle position data for seed ``XXXXXX``
        * ``XXXXXX-prof`` raw density profile data for seed ``XXXXXX``
        * ``XXXXXX-param`` parameter info for seed ``XXXXXX``
    * ``0.125/`` discretization $\alpha=0.125$ data
        * same content as ``0.0/`` but for $\alpha=0.125$
    * etc.