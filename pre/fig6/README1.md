# Revisiting the Ratchet Principle companion, fig 6 #

This directory contains code and data necessary for generating Fig. 6.

Thhe figure is generated using the script ``fig6.py``, which uses the parsed data in the ``data/`` folder.

* ``fig6.py`` Code to analyze data and generate figure
* ``data/``
    * ``N12/``
        * ``J_avg`` Averaged current for all simulations with $N=12$. Leftmost column: $N$. Second column: $N*\varepsilon$, where $\varepsilon=\frac{2}{3} \text{amp}$. Third column: Average current divided by system size (averaged over seeds and particles). Fourth column: Standard deviation of current divided by $L$ (taken across the ensemble of seeds after each one is particle-averaged).
        * ``amp-0.375-avg_disp`` particle displacement vs. time averaged over all seeds with interaction strength ``amp=0.375``. Leftmost column is time, middle column is average displacement divided by system size $L$, and rightmost column is the standard error of the avg displacement/$L$ over the ensemble of seeds.
        * ``amp-0.375-avg_profile_rho`` average density profile for all simulations with interaction strength ``amp=0.375``. Leftmost column is $x$, and rightmost column is the average density. Needs to be normalized before plotting.
        * ``amp-0.375-avg_profile_m`` average magnetization profile for all simulations with interaction strength ``amp=0.375``. Leftmost column is $x$, and rightmost column is the average magnetization. Needs to be normalized before plotting.
        * ``amp_0.375-1348-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1348``. The first line is the command used in the simulation.
        * ``amp_0.375-1705-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1705``. The first line is the command used in the simulation.
        * etc.
    * ``N25/``
        * etc.
    * etc.