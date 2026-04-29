# Revisiting the Ratchet Principle (companion), fig 1 #

This directory contains code and data necessary for generating Fig. 1.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash script ``C_code/loop.sh``. It generates data that it saves in the folder ``data``. A copy of the data is already there.

The data is parsed using the script ``data/get_data.py`` and left in the folders ``data/ABPs/``, etc. A copy of the parsed data is already there.

Finally, the figure is generated using the script ``fig1.py``, which uses the parsed data in the ``data/`` folder.

* ``fig1.py`` Code to create figure
* ``C_code/``
    * ``loop.sh`` Bash script to run ABP, PBP, AOUP, and UPBP simulations in C. Runs 8 different type of simulation (interaction and non-interacting of each system), 15 seeds for each
    * ``util/``
        * ``inhom_fluct.c`` functions to define cubic inhomogeneous fluctuation landscape
        * ``PFs.c`` functions for calculating pairwise forces
        * ``spatial_hashing.c`` functions for looping through particle force calculations using spatial hashing algorithm
        * ``util.c`` miscellaneous functions
        * ``mt64.h`` For the random-number generator
        * ``mt19937-64.c`` For the random-number generator
    * ``ABPs/``
        * ``ABP-RTP-2D.c`` Main simulation code for ABPs
        * ``ABP-RTP-2D-functions.c`` Functions for the ABP simulations
        * ``options.h`` Pre-compiler options for the ABP simulations
        * ``abp-cub-disp-prof-J`` Compiled interacting ABP simulator
        * ``abp-ni-disp-prof-J`` Compiled non-interacting ABP simulator
    * ``AOUPs/``
        * ``AOUPs-2D.c`` Main simulation code for AOUPs
        * ``AOUPs-2D-functions.c`` Functions for the AOUP simulations
        * ``options.h`` Pre-compiler options for the AOUP simulations
        * ``aoup-cub-disp-prof-J`` Compiled interacting AOUP simulator
        * ``aoup-ni-disp-prof-J`` Compiled non-interacting AOUP simulator
    * ``PBPs/``
        * ``OD-PBPs-2D.c`` Main simulation code for overdamped PBPs
        * ``OD-PBPs-2D-functions.c`` Functions for the PBP simulations
        * ``options.h`` Pre-compiler options for the PBP simulations
        * ``pbp-cub-disp-prof-J`` Compiled interacting PBP simulator
        * ``pbp-ni-disp-prof-J`` Compiled non-interacting PBP simulator
    * ``UPBPs/``
        * ``UD-PBPs-2D.c`` Main simulation code for underdamped PBPs (UPBPs)
        * ``UD-PBPs-2D-functions.c`` Functions for the UPBP simulations
        * ``options.h`` Pre-compiler options for the UPBP simulations
        * ``upbp-cub-disp-prof-J`` Compiled interacting UPBP simulator
        * ``upbp-ni-disp-prof-J`` Compiled non-interacting UPBP simulator
* ``data/``
    * ``get_data.py`` Gets the results from the C simulations and outputs in the format needed for the ``fig1.py`` script
    * ``ABPs/``
        * ``int/``
            * ``run.sh`` script that runs simulation of interacting ABPs with provided parameters
            * ``tf_8000000-dt_0.000250-N_120-Lx_30-Ly_10-Nmode_7-k_10/``
                * ``prof_avg`` averaged density profile over all seeds in this folder
                * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
                * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
                * ``Jxprof_avg`` averaged x-current density profile over all seeds in this folder
                * ``Jyprof_avg`` averaged y-current density profile over all seeds in this folder
                * ``XXXXX-param`` parameter info for seed ``XXXXX``
                * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
                * ``XXXXX-prof`` density profile data for seed ``XXXXX``
                * ``XXXXX-Jprof`` current density profile data for seed ``XXXXX``
                * ``XXXXX-pos`` particle position data for seed ``XXXXX``
        * ``nonint/``
            * same contents as ``int/`` but for non-interacting simulations
    * ``AOUPs/``
        * same contents as ``ABPs/`` but for AOUPs
    * ``UPBPs/``
        * same contents as ``ABPs/`` but for UPBPs
    * ``PBPs/``
        * same contents as ``ABPs/`` but for PBPs
