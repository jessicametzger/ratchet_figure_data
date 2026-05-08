#!/bin/bash

vars=("3958"
"3115"
"3839"
"3844"
"3619"
"3926"
"3919"
"3548"
)
for var_1 in "${vars[@]}"; do
  nohup ./run.sh $var_1 &
done

wait

