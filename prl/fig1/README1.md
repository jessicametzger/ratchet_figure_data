# Revisiting the Ratchet Principle, fig 1 #

This directory contains code and data necessary for generating Fig. 1.

The figure is generated using the script ``fig1.py``, which uses the parsed data in the ``data/`` folder.

* ``fig1.py`` Code to analyze data
* ``data/``
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
