# Revisiting the Ratchet Principle (companion), fig 10 #

This directory contains code and data necessary for generating Fig. 10.

The figure is generated using the script ``fig10.py``, which uses the data in the ``data/asymm_T-step_V/``, ``data/symm_T-symm_V/``, and ``data/step_T-symm_V/`` folders.

* ``fig10.py`` Code to create figure
* ``data/``
    * ``asymm_T-step_V/``
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