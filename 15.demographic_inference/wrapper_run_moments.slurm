#!/bin/bash

#SBATCH --mem 2G
#SBATCH --time 00:30:00
#SBATCH --partition standard
#SBATCH --account berglandlab_standard
#SBATCH -c 1
#SBATCH -e slurm/%a.err
#SBATCH --output slurm/%a.out

options_file=options/options_run_moments_Australia.tsv
script_file=run_moments.py

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
