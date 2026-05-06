#!/bin/bash

# ./aoups-cub-disp-prof-sig file Lx Ly N dt tf seed tau Nstep_D Dxs Ds a k StoreInterPos StoreInterDisp Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf
../../C_code/aoups-cub-disp-prof-sig $1 60 15 700 0.00025 20000.0 $1 1 4 0,10,20,50 20,1,1,20 1 50 20000 2000 600 1 20000 20000 1
