#!/bin/bash

U1=$(echo "scale=8; $1 * 19" | bc)


# ./pbp-Tofx-Uofx file Lx Ly N dt tf seed Nstep_T T_xs Ts Nstep_U U_xs Us a k StoreInterPos StoreInterDisp Nbinx Nbiny NstepProf StoreInterProf UpdateInterProf
../C_code/pbp-Tofx-Uofx $1/$2 20 5 200 0.00025 20000 $2 3 0,2.5,20 20,1,20 3 0,2.5,20 0,$U1,0 1 50 20000 2000 200 1 20000 20000 1