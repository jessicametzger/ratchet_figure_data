# Revisiting the Ratchet Principle (companion), fig 7 #

This directory contains code and data necessary for generating Fig. 7.

The figure is generated using the script ``fig7.py``, which uses the parsed data in the ``data/int/`` and ``data/nonint/`` folders.

* ``fig7.py`` Code to create figure
* ``data/``
    * ``int/``
        * ``prof_avg`` averaged density profile over all seeds in this folder
        * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
        * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
        * ``sigm_avg`` averaged active stress tensor over all seeds in this folder. Each row is a different x point (evenly spaced from 0 to Lx=60). Column 1 is $\rho D/\tau$, column 2 is $\langle v_x^2\rangle$, column 3 is $\langle v_x v_y\rangle$, column 4 is $\langle v_y^2\rangle$, column 5 is $\langle v_x F_x\rangle$, column 6 is $\langle v_x F_y\rangle$, column 7 is $\langle v_y F_x\rangle$, and column 8 is $\langle v_y F_y\rangle$.
        * ``sigIK_avg`` averaged Irving-Kirkwood interaction stress over all seeds in this folder. Each row is a different x point (evenly spaced from 0 to Lx=60). Column 1 is $ \sigma^{\rm IK}_{xx}$, column 2 is $ \sigma^{\rm IK}_{xy}$, and column 3 is $\sigma^{\rm IK}_{yy}$.
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-sigmaAprof`` thermal/active stress density data for seed ``XXXXX``
        * ``XXXXX-sigmaIKprof`` interaction (Irving-Kirkwood) stress density data for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
    * ``nonint/``
        * ``prof_avg`` averaged density profile over all seeds in this folder
        * ``disp_avg`` averaged particle displacement over time over all seeds in this folder
        * ``disp_err`` standard error of the mean particle displacement over time over all seeds in this folder
        * ``XXXXX-param`` parameter info for seed ``XXXXX``
        * ``XXXXX-disp`` particle displacement data for seed ``XXXXX``
        * ``XXXXX-prof`` density profile data for seed ``XXXXX``
        * ``XXXXX-sigmaAprof`` thernal/active stress density data for seed ``XXXXX``
        * ``XXXXX-pos`` particle position data for seed ``XXXXX``
