# Revisiting the Ratchet Principle companion, fig 5 #

This directory contains code and data necessary for generating Fig. 5.

The data is saved in the ``data/`` directory. There are 10 subdirectories in here, each of the form ``data/N_500-amp_0.XX-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12/`` corresponding to a different interaction strength ``0.XX``. To generate the data for that interaction strength, run the script ``data/N_500-amp_0.XX-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12/loop.sh``. It executes the C code in the directory ``C_code/``. A copy of the resulting data is stored there.

Then, the figure is generated using the script ``fig5.py``, which uses the data in the ``data/`` folder.

* ``fig5.py`` Code to analyze data
* ``C_code/``
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``Ratchet-PFAPs-2D.c`` Main code for ABP/RTP simulations
    * ``Ratchet-PFAPs-2D-functions.c`` Functions for ABP/RTP simulations
    * ``spatial_hashing.c`` Methods for spatial hashing algorithm used to treat interactions
    * ``ratchet-inv_h-abp-cubic-force-sigma`` compiled C code (for ABPs)
    * ``ratchet-inv_h-rtp-cubic-force-sigma`` compiled C code (for RTPs)
    * ``options.h`` pre-compiler options for C code (e.g. ABP vs. RTP, interacting vs. non-interacting)
* ``data/``
    * ``N_500-amp_0.1-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12/`` data for simulations with interaction strength 0.1
        * ``loop.sh`` bash script to run simulations in this folder (for all seeds, both ABPs and RTPs)
        * ``run.sh`` bash script to execute simulation with provided seed, for provided system (ABP or RTP). Called by ``loop.sh``.
        * ``abp_1021/`` data for ABP simulations with seed 1021
            * ``res-param`` parameter info for this run
            * ``res-disp`` particle net displacement data (columns: time, particleID, x displacement, y displacement)
            * ``res-pos`` particle position data (columns: time, particleID, x, y, theta)
            * ``res-profile_rho`` density profile, given by net particle counts in a list of boxes (columns: time, x, y, particle count)
        * ``abp_1117/`` data for ABP simulations with seed 1117
        * etc.
        * ``rtp_1093/`` data for RTP simulations with seed 1093
        * ``rtp_1183/`` data for RTP simulations with seed 1183
        * etc.
    * ``N_500-amp_0.2-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12/`` data for simulations with interaction strength 0.2
    * etc.