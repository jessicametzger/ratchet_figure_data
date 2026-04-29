#!/bin/bash

# ./upbps-ni-cub-disp-prof-J-sig file Lx Ly N dt tf seed gamma Nstep_T Txs Ts StoreInterPos StoreInterDisp Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf
../../C_code/upbps-ni-cub-disp-prof-J-sig $1 60 15 700 0.00025 200.0 $1 1 4 0,10,20,50 20,1,1,20 200 20 600 1 200 200 1
