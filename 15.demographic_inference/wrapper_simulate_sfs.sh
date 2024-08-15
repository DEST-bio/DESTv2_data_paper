#!/bin/bash

cli_params=`sed "${1}q;d" options/options_simulate_sfs.tsv`

python simulate_sfs.py $cli_params
