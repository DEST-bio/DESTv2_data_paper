#!/bin/bash

#SBATCH --mem 10G
#SBATCH --time 6:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab_standard
#SBATCH -c 1
#SBATCH -e slurm/%a.err
#SBATCH --output slurm/%a.out

options_file=options/options_construct_jackknife_sfs.txt
script_file=construct_jackknife_sfs.py

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
