#!/bin/bash

for i in {0..39}; do
  python run_moments_on_simulated_sfss.py $1 $2 $i &
done

wait
