#!/bin/bash

vars=("2109"
"2051"
"2420"
"2410"
"2647"
"2111"
"2959"
"2105"
)
for var_1 in "${vars[@]}"; do
  nohup ./run.sh $var_1 &
done

wait

