#!/bin/bash

vars=("2746"
"2852"
"2745"
"2206"
"2437"
"2203"
"2251"
"2396"
)
for var_1 in "${vars[@]}"; do
  nohup ./run.sh $var_1 &
done

wait

