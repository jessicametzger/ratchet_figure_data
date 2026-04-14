#!/bin/bash

# run simulation code with the provided seed
./abp-cub-sigA-sigIK-FA-Fint ../data/$1 40 15 400 0.0005 700000 1 4 0.0,4.0,14.0,30.0 1.0,30.0,30.0,1.0 1 50 $1 700000 400 1 699900 700000 1 14000
