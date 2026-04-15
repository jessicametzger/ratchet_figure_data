# Revisiting the Ratchet Principle, fig 2 #

This directory contains code and data necessary for generating Fig. 2.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The simulation is run by executing the following command from this directory:

``$ ./C_code/pbp-sigIK-epr ./data/1337 50 15 400 0.0005 100000 1337 4 0,8,17,41 10,1,1,10 1 50 100000 500 150 20000000 10000 0.0005``

It saves the data in the files ``data/1337-...``. Due to space constraints, a copy of the raw data is not stored here.

Then, the data is parsed and averaged using the script ``avg_data.py``. The averaged data files are stored in the ``data/`` folder. A copy of the parsed and averaged data is included here.

Finally, the figure is generated using the script ``fig2.py``, which uses the parsed data in the ``data/`` folder.

* ``fig2.py`` Code to analyze data
* ``C_code/``
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``OD-PBPs-2D.c`` Main code for PBP simulations
    * ``OD-PBPs-2D-functions.c`` Functions for PBP simulations
    * ``spatial_hashing.c`` Methods for spatial hashing algorithm used to treat interactions
    * ``inhom_fluct.c`` Methods for computing activity landscape
    * ``PFs.c`` Methods for computing pairwise forces and interaction stresses
    * ``pbp-sigIK-epr`` compiled C code
* ``data/``
    * ``avg_data.py`` Code to parse and average data
    * ``eprprof_avg`` The averaged EPR heatmap. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$. Each column is a different $y$ value, equally spaced from 0 to $L_y=15$.
    * ``eprseries_avg`` Time series of the cumulative EPR. Each row is a different time point, equally spaced from 0 to $t_f=10^5$.
    * ``sigprof_avg`` The average stress profile. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
        * Column 1 contains the $x$ values 
        * Column 2 contains the average of $\sigma^{\rm IK}_{xx}$ 
        * Column 3 contains the average of $\sigma^{\rm IK}_{xy}=\sigma^{\rm IK}_{yx}$ 
        * Column 4 contains the average of $\sigma^{\rm IK}_{yy}$
        * Column 5 contains the average of $\sum_i T(\mathbf{r}_i)$
    * ``1337-param`` File containing the parameter info for a simulation