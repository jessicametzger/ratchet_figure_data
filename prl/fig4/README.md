# Revisiting the Ratchet Principle, fig 4 #

This directory contains code and data necessary for generating Fig. 4.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The simulation is run using the bash script ``C_code/loop.sh``, which repeatedly calls the script ``run.sh`` for different alpha (tumbling rate) and epsilon (potential strength).

It saves the data in the folders ``data/alphaX.00/epsilonX.XXX/res-...``. A copy of the data is included there.

Then, the data is parsed and averaged using the script ``data/avg_data.py``. The averaged data files are stored in the ``data/`` folder. A copy of the parsed and averaged data is included here.

Finally, the figure is generated using the script ``fig4.py``, which uses the parsed data in the ``data/`` folder.

* ``fig4.py`` Code to analyze data
* ``C_code/``
    * ``loop.sh`` script which repeatedly calls the script ``run.sh`` with different parameters to run the RTP simulation
    * ``run.sh`` script which, given 2 arguments alpha and epsilon, runs the RTP simulation with those values of alpha and epsilon
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``Ratchet-potential-activity-simple-1.c`` Main code for RTP simulations
    * ``Ratchet-potential-activity-simple-1-functions.c`` Functions for RTP simulations
    * ``simple-1`` compiled C code
* ``data/``
    * ``avg_data.py`` Code to parse and average data
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