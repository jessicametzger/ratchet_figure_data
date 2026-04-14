# Revisiting the Ratchet Principle, fig 3 #

This directory contains code and data necessary for generating Fig. 3.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The compiled code is looped by the bash scripts ``loop.sh`` and ``run.sh``. It generates data that it saves in the files ``data/[seed]-...``. 

Then, the data is parsed and averaged using the script ``avg_data.py``. The averaged data files are stored in the ``data/`` folder.

Finally, the figure is generated using the script ``fig3.py``, which uses the parsed data in the ``data/`` folder.

* ``fig3.py`` Code to analyze data
* ``C_code/``
    * ``loop.sh`` Bash script to run ABP simulations in C. Runs 128 simulations, each with a different random number generator seed, by calling ``run.sh``. Places the results in the ``data/`` folder.
    * ``run.sh`` Script to run ABP simulation in C
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``ABP-v_step-2D.c`` Main code for ABP simulations
    * ``ABP-v_step-2D-functions.c`` Functions for ABP simulations
    * ``spatial_hashing.c`` Methods for spatial hashing algorithm used to treat interactions
    * ``inhom_fluct.c`` Methods for computing activity landscape
    * ``PFs.c`` Methods for computing pairwise forces and interaction stresses
    * ``abp-v_step-2d`` compiled C code
* ``data/``
    * ``avg_data.py`` Code to parse and average data
    * ``disp_avg`` The average particle x-displacement at equally-spaced times up to tf=$7 \times 10^5$. Averaged over seeds, summed over particles.
    * ``FAprof_avg`` The average active force profile. Each row is a different $x$ value, equally spaced from 0 to $L_x=40.$
        * Column 1 contains the average of $\sum_i [\dot{x}_i \partial_x v(\mathbf{r}_i) + \dot{y}_i \partial_y v(\mathbf{r}_i)] \cos\theta_i$
        * Column 2 contains the average of $\sum_i [\dot{x}_i \partial_x v(\mathbf{r}_i) + \dot{y}_i \partial_y v(\mathbf{r}_i)] \sin\theta_i$
    * ``Fintprof_avg`` The average interaction force profile. Each row is a different $x$ value, equally spaced from 0 to $L_x=40.$
        * Column 1 contains the average of $\sum_i F^x(\mathbf{r}_i)$
        * Column 2 contains the average of $\sum_i F^y(\mathbf{r}_i)$
    * ``sigAprof_avg`` The average active stress profile.  Each row is a different $x$ value, equally spaced from 0 to $L_x=40.$
        * Column 1 contains the average of $\sum_i \dot{x}_i v(\mathbf{r}_i) \tau \cos\theta_i$
        * Column 2 contains the average of $\sum_i \dot{y}_i v(\mathbf{r}_i) \tau \cos\theta_i$
        * Column 3 contains the average of $\sum_i \dot{x}_i v(\mathbf{r}_i) \tau \sin\theta_i$
        * Column 4 contains the average of $\sum_i \dot{y}_i v(\mathbf{r}_i) \tau \sin\theta_i$
    * ``sigIKprof_avg`` The average interaction (Irving-Kirkwood) stress. Each row is a different $x$ value, equally spaced from 0 to $L_x=40.$
        * Column 1 contains the average of $\sigma^{\rm IK}_{xx}$ 
        * Column 2 contains the average of $\sigma^{\rm IK}_{xy}=\sigma^{\rm IK}_{yx}$ 
        * Column 3 contains the average of $\sigma^{\rm IK}_{yy}$
    * ``xxxxx-param`` File containing the parameter info for simulation with seed ``xxxxx``.
    * ``xxxxx-disp`` Raw displacement data for simulation with seed ``xxxxx``. Even rows are the $x$ displacement, and odd rows are the $y$ displacement. The leftmost column is the time at measurement. The $(i+1)$th column is the displacement of the $i$th particle.
    * ``xxxxx-FAprof`` Raw active force profile data for simulation with seed ``xxxxx``. Leftmost column is the time of recording. Right column is the data. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
    * ``xxxxx-Fintprof`` Raw interaction force profile data for simulation with seed ``xxxxx``. Left column is the time of recording. Right column is the data. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
    * ``xxxxx-pos`` File containing list of positions of particles in simulation with seed ``xxxxx`` at end of simulation. Row 1 columns 2-401 contains the $x$ positions; row 2 columns 2-401 contains the $y$ positions.
    * ``xxxxx-prof`` File containing the density profile of particles in simulation with seed ``xxxxx``. Leftmost column is the time of recording. Right column is the data. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
    * ``xxxxx-sigmaAprof`` File containing the active stress profile of particles in simulation with seed ``xxxxx``. Leftmost column is the time of recording. Right column is the data. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
    * ``xxxxx-sigmaIKprof`` File containing the interaction (Irving-Kirkwood) stress of simulation with seed ``xxxxx``. Leftmost column is the time of recording. Rows 1-40 contain the data for $\sigma^{\rm IK}_{xx}$. Rows 41-80 contain the data for $\sigma^{\rm IK}_{xy}=\sigma^{\rm IK}_{yx}$ . Rows 81-120 contain the data for $\sigma^{\rm IK}_{yy}$. Within each section (rows 1-40, 41-80, and 81-120), the different rows correspond to different $x$-values, evenly spaced from 0 to $L_x=40$. Each column from columns 2-16 corresponds to a different $y$ value, evenly spaced from 0 to $L_y=15$. Note that this data is on a coarser grid, since it is measured differently than the other profiles.