# Revisiting the Ratchet Principle (companion), fig 1 #

This directory contains code and data necessary for generating Fig. 1.

The figure is generated using the script ``fig1.py``, which uses the parsed data in the ``data/`` folder.

* ``fig1.py`` Code to create figure
* ``data/``
    * ``ABPs/``
        * ``int/``
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
