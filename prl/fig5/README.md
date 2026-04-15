# Revisiting the Ratchet Principle, fig 5 #

This directory contains code and data necessary for generating Fig. 5.

The data is generated through simulation of the C code contained in the directory ``C_code/``. The simulations are run using the bash scripts ``C_code/N[...]_run_sims.sh``, each of which which repeatedly calls the script ``ratchet-pfaps`` for that $N$ with different interaction strength and random number seed.

It saves the data for ``N=X`` with interaction strength ``amp=0.yyy`` and random number seed ``zzzzz`` in files prefixed by ``data/NX/amp_0.yyy-zzzzz-``. Due to space constraints, we have only included the ``-param`` file for each simulation, containing a list of parameters used in that simulation.

Then, the data is parsed and averaged using the script ``data/avg_data.py``. The averaged current for ``N=X`` is stored in the file ``data/NX/J_avg``. A copy of the parsed and averaged data is included here for each $N$.

Finally, the figure is generated using the script ``fig5.py``, which uses the parsed data in the ``data/`` folder.

* ``fig5.py`` Code to analyze data and generate figure
* ``mean_field_schematic_thin.pdf`` pre-made schematic that is included as panel (a) of the figure by ``fig5.py``
* ``C_code/``
    * ``N12_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=12$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``N25_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=25$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``N50_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=50$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``N100_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=100$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``N200_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=200$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``N400_run_sims.sh`` script which repeatedly runs ``C_code/ratchet-pfaps`` to simulate $N=400$ RTPs with different interaction strengths ``amp``, and different random number seeds.
    * ``mt64.h`` For the random-number generator
    * ``mt19937-64.c`` For the random-number generator
    * ``Ratchet-PFAPs.c`` Main code for RTP simulations
    * ``Ratchet-PFAPs-functions.c`` Functions for RTP simulations
    * ``ratchet-pfaps`` compiled C code
* ``data/``
    * ``avg_data.py`` Code to parse and average current data
    * ``N12/``
        * ``J_avg`` Averaged current for all simulations with $N=12$. Leftmost column: $N$. Second column: $N*\varepsilon$, where $\varepsilon=\frac{2}{3} \text{amp}$. Third column: Average current divided by system size (averaged over seeds and particles). Fourth column: Standard deviation of current divided by $L$ (taken across the ensemble of seeds after each one is particle-averaged).
        * ``amp_0.375-1348-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1348``. The first line is the command used in the simulation.
        * ``amp_0.375-1705-param`` File containing parameter info for the $N=12$ simulation with interaction strength ``amp=0.375`` and random number seed ``1705``. The first line is the command used in the simulation.
        * etc.
    * ``N25/``
        * etc.
    * etc.