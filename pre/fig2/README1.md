# Revisiting the Ratchet Principle (companion), fig 2 #

This directory contains code and data necessary for generating Fig. 2.

The figure is generated using the script ``fig2.py``, which uses the parsed data in the ``data/int/`` and ``data/nonint/`` folders.

* ``fig2.py`` Code to create figure
* ``data/``
    * ``int/``
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
        * ``prof_avg`` averaged density profile over all seeds in this folder
        * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
        * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-sigmaprof`` thernal stress density data for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
