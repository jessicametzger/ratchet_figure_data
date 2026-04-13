# Revisiting the Ratchet Principle, fig 1 #

This directory contains code and data necessary for generating Fig. 1.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash script ``run_sims.sh``. It generates data that it saves in the files ``ABP-int``, ``ABP-ni``, ``PBP-int``, and ``PBP-ni``. A copy of the data is already there.

The data is parsed using the script ``data/get_data.py`` and left in the folders ``data/ABPs/`` and ``data/PBPs/``. A copy of the parsed data is already there.

Finally, the figure is generated using the script ``fig1.py``, which uses the parsed data in the ``data/`` folder.

* ``fig1.py`` Code to analyze data
* ``C_code/``
    * ``run_sims.sh`` Bash script to run ABP and PBP simulations in C. Runs 4 simulations: ABP interacting, ABP non-interacting, PBP interacting, and PBP non-interacting
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``ABP/``
        * ``Ratchet-PFAPs-2D.c`` Main simulation code for ABPs
        * ``Ratchet-PFAPs-2D-functions.c`` Functions for the ABP simulations
        * ``options.h`` Pre-compiler options for the ABP simulations
        * ``abp-int`` Compiled interacting ABP simulator
        * ``abp-ni`` Compiled non-interacting ABP simulator
    * ``PBP/``
        * ``PBPs-2D.c`` Main simulation code for PBPs
        * ``PBPs-2D-functions.c`` Functions for the PBP simulations
        * ``options.h`` Pre-compiler options for the PBP simulations
        * ``pbp-int`` Compiled interacting PBP simulator
        * ``pbp-ni`` Compiled non-interacting PBP simulator
* ``data/``
    * ``get_data.py`` Gets the results from the C simulations and outputs in the format needed for the ``fig1.py`` script
    * ``ABPs/``
        * ``active-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025-disp_avg`` Particle-averaged net displacement (active, interacting)
        * ``active-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025-prof_avg`` Average density profile (active, interacting)
        * ``active-ni-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025-disp_avg`` Particle-averaged net displacement (active, non-interacting)
        * ``active-ni-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025-prof_avg`` Average density profile (active, non-interacting)
        * ``active-param``
    * ``PBPs/``
        * ``passive-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025-disp_avg`` Particle-averaged net displacement (passive, interacting)
        * ``passive-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025-prof_avg`` Average density profile (passive, interacting)
        * ``passive-ni-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025-disp_avg`` Particle-averaged net displacement (passive, non-interacting)
        * ``passive-ni-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025-prof_avg`` Average density profile (passive, non-interacting)
        * ``passive-param``
