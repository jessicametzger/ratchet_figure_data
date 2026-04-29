# Revisiting the Ratchet Principle, fig 4 #

This directory contains code and data necessary for generating Fig. 4.

The figure is generated using the script ``fig4.py``, which uses the parsed data in the ``data/`` folder.

* ``fig4.py`` Code to analyze data
* ``data/``
    * ``current_X.00`` The average current observed in the simulations with alpha=X.00, for different epsilon.
        * Column 1 contains the value of epsilon
        * Column 2 contains the average current seen in simulations
    * ``current_th_X.00`` The current predicted by our theory for alpha=X.00, for different epsilon.
        * Column 1 contains the value of epsilon
        * Column 2 contains the predicted current
    * ``alphaX.00/`` folder containing simulation data with tumbling rate alpha=X.00
        * ``epsilonX.XXX/`` folder containing a simulation with epsilon=X.XXX
            * ``res-disp`` net displacement (rightmost column) for each particle (second column) at different times (leftmost column)
            * ``res-param`` parameter info for the simulation
            * The rest of the files in this folder are empty (they are not used in the analysis).