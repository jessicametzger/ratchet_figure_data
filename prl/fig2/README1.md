# Revisiting the Ratchet Principle, fig 2 #

This directory contains code and data necessary for generating Fig. 2.

The figure is generated using the script ``fig2.py``, which uses the parsed data in the ``data/`` folder.

* ``fig2.py`` Code to analyze data
* ``data/``
    * ``eprprof_avg`` The averaged EPR heatmap. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$. Each column is a different $y$ value, equally spaced from 0 to $L_y=15$.
    * ``eprseries_avg`` Time series of the cumulative EPR. Each row is a different time point, equally spaced from 0 to $t_f=10^5$.
    * ``sigprof_avg`` The average stress profile. Each row is a different $x$ value, equally spaced from 0 to $L_x=40$.
        * Column 1 contains the $x$ values 
        * Column 2 contains the average of $\sigma^{\rm IK}_{xx}$ 
        * Column 3 contains the average of $\sigma^{\rm IK}_{xy}=\sigma^{\rm IK}_{yx}$ 
        * Column 4 contains the average of $\sigma^{\rm IK}_{yy}$
        * Column 5 contains the average of $\sum_i T(\mathbf{r}_i)$
    * ``1337-param`` File containing the parameter info for a simulation