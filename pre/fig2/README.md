# Revisiting the Ratchet Principle (companion), fig 2 #

This directory contains code and data necessary for generating Fig. 2.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash scripts ``data/int/loop.sh`` (interacting) and ``data/nonint/loop.sh`` (non-interacting). It generates data that it saves in the folders ``data/int/`` and ``data/nonint/``. A copy of the data is already there.

The data is parsed using the script ``data/get_data.py`` and left in the folders ``data/int/`` and ``data/nonint/``. A copy of the parsed data is already there.

Finally, the figure is generated using the script ``fig2.py``, which uses the parsed data in the ``data/int/`` and ``data/nonint/`` folders.

* ``fig2.py`` Code to create figure
* ``C_code/``
    * ``util/``
        * ``inhom_fluct.c`` functions to define cubic inhomogeneous fluctuation landscape
        * ``PFs.c`` functions for calculating pairwise forces
        * ``spatial_hashing.c`` functions for looping through particle force calculations using spatial hashing algorithm
        * ``util.c`` miscellaneous functions
        * ``mt64.h`` For the random-number generator
        * ``mt19937-64.c`` For the random-number generator
    * ``UD-PBPs-2D.c`` Main simulation code for UD PBPs
    * ``UD-PBPs-2D-functions.c`` Functions for the UD PBP simulations
    * ``options.h`` Pre-compiler options for the UD PBP simulations
    * ``abp-cub-disp-prof-J`` Compiled interacting UD PBP simulator
    * ``abp-ni-disp-prof-J`` Compiled non-interacting UD PBP simulator
* ``data/``
    * ``get_data.py`` Gets the results from the C simulations and outputs in the format needed for the ``fig2.py`` script
    * ``int/``
        * ``run.sh`` script that runs simulation of interacting UD PBPs with provided parameters
        * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
        * ``prof_avg`` averaged density profile over all seeds in this folder
        * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
        * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
        * ``sigm_avg`` averaged thermal stress tensor over all seeds in this folder. Each row is a different x point (evenly spaced from 0 to Lx=60). Column 1 is $\langle v_x^2\rangle$, column 2 is $\langle v_x v_y\rangle$, and column 3 is $\langle v_y^2\rangle$.
        * ``sigIK_avg`` averaged Irving-Kirkwood interaction stress over all seeds in this folder. Each row is a different x point (evenly spaced from 0 to Lx=60). Column 1 is $ \sigma^{\rm IK}_{xx}$, column 2 is $ \sigma^{\rm IK}_{xy}$, and column 3 is $\sigma^{\rm IK}_{yy}$.
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-sigmaprof`` thernal stress density data for seed ``XXXXX``
        * ``XXXXX-sigmaIKprof`` interaction (Irving-Kirkwood) stress density data for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
    * ``nonint/``
        * ``run.sh`` script that runs simulation of non-interacting UD PBPs with provided parameters
        * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
        * ``prof_avg`` averaged density profile over all seeds in this folder
        * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
        * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-sigmaprof`` thernal stress density data for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
