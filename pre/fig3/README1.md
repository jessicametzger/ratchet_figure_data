# Revisiting the Ratchet Principle (companion), fig 3 #

This directory contains code and data necessary for generating Fig. 3.

The figure is generated using the script ``fig3.py``, which uses the parsed data in the ``data/`` folder.

* ``fig3.py`` Code to create figure
* ``data/``
    * ``0.0/`` discretization $\alpha=0$ data
        * ``disp_avg`` averaged particle displacements for evenly-spaced time intervals
        * ``disp_err`` standard error of the mean average particle displacement over all seeds, for evenly-spaced time intervals
        * ``prof_avg`` average particle density for evenly spaced $x$ from 0 to Lx
        * ``XXXXXX-disp`` raw particle displacement data for seed ``XXXXXX``
        * ``XXXXXX-pos`` raw particle position data for seed ``XXXXXX``
        * ``XXXXXX-prof`` raw density profile data for seed ``XXXXXX``
        * ``XXXXXX-param`` parameter info for seed ``XXXXXX``
    * ``0.125/`` discretization $\alpha=0.125$ data
        * same content as ``0.0/`` but for $\alpha=0.125$
    * etc.