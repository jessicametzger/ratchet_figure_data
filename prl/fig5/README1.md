# Revisiting the Ratchet Principle, fig 5 #

This directory contains code and data necessary for generating Fig. 5.

Thhe figure is generated using the script ``fig5.py``, which uses the parsed data in the ``data/`` folder.

* ``fig5.py`` Code to analyze data and generate figure
* ``mean_field_schematic_thin.pdf`` pre-made schematic that is included as panel (a) of the figure by ``fig5.py``
* ``data/``
    * ``N12/``
        * ``J_avg`` Averaged current for all simulations with $N=12$. Leftmost column: $N$. Second column: $N*\varepsilon$, where $\varepsilon=\frac{2}{3} \text{amp}$. Third column: Average current divided by system size (averaged over seeds and particles). Fourth column: Standard deviation of current divided by $L$ (taken across the ensemble of seeds after each one is particle-averaged).
        * ``amp_0.375-1348-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1348``. The first line is the command used in the simulation.
        * ``amp_0.375-1705-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1705``. The first line is the command used in the simulation.
        * etc.
    * ``N25/``
        * etc.
    * etc.