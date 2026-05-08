#!/bin/bash

vars=("3178"
"3536"
"3165"
"3531"
"3461"
"3926"
"3358"
"3056"
)
for var_1 in "${vars[@]}"; do
  nohup ./run.sh $var_1 &
done

wait

