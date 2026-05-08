#!/bin/bash

vars=("3732"
"3943"
"3361"
"3811"
"3155"
"3330"
"3841"
"3289"
)
for var_1 in "${vars[@]}"; do
  nohup ./run.sh $var_1 &
done

wait

