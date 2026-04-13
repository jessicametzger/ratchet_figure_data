#!/bin/bash

# run simulations in parallel

# usage: ./abp-int file Lx Ly N dt tf tau Nstep_v v_x1,v_x2,... v1,v2,... a k seed StoreInterPos Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf StoreInterDisp
ABP/abp-int ABP-int 80 50 800 0.00025 1000000.0 1 5 0,40,42.5,62.5,75 20,20,1,1,20 1 50 9826 1000000 800 1 500000 1000000 2 10000 &

# usage: ./abp-ni  file Lx Ly N dt tf tau Nstep_v v_x1,v_x2,... v1,v2,... seed StoreInterPos Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf StoreInterDisp
ABP/abp-ni  ABP-ni  80 50 800 0.00025 1000000.0 1 5 0,40,42.5,62.5,75 20,20,1,1,20 2269 1000000 800 1 500000 1000000 2 10000 &

# usage: ./pbp-int file Lx Ly N dt tf seed Nstep_T T_xs Ts a k StoreInterPos StoreInterDisp Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf
PBP/pbp-int PBP-int 80 50 800 0.00025 1000000.0 8761 5 0,40,42.5,62.5,75 20,20,1,1,20 1 50 1000000 10000 800 1 500000 1000000 2 &

# usage: ./pbp-ni  file Lx Ly N dt tf seed Nstep_T T_xs Ts StoreInterPos StoreInterDisp Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf
PBP/pbp-ni  PBP-ni  80 50 800 0.00025 1000000.0 4976 5 0,40,42.5,62.5,75 20,20,1,1,20 1000000 10000 800 1 500000 1000000 2
