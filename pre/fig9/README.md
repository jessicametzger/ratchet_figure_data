# Revisiting the Ratchet Principle (companion), fig 9 #

This directory contains code and data necessary for generating Fig. 9.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The simulation is run using the bash script ``data/loop.sh``, which repeatedly calls the script ``data/run.sh`` for different random number seed. It saves the data in ``data/``. A copy of the data is included there.

The figure is generated using the script ``fig9.py``, which uses the data in the ``data/`` folder.

* ``fig9.py`` Code to analyze data
* ``C_code/``
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``options.h`` Pre-compiler options for C code
    * ``Ratchet-potential-activity-trig.c`` Main code for RTP simulations
    * ``Ratchet-potential-activity-trig-functions.c`` Functions for RTP simulations
    * ``ratchet`` compiled C code
* ``data/``
    * ``run.sh`` script that runs ``C_code/ratchet`` with provided random number seed
    * ``loop.sh`` script that repeatedly calls ``run.sh`` for different random number seeds
    * ``XXXXX-param`` parameter info for seed ``XXXXX``
    * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
    * ``XXXXX-profile_rho`` density profile data for seed ``XXXXX``
    * ``XXXXX-profile_m`` magnetization profile data for seed ``XXXXX``