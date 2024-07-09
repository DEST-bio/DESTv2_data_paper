#!/bin/bash

#SBATCH --mem 3G
#SBATCH --time 01:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab_standard
#SBATCH -c 1
#SBATCH -e slurm/%a.err
#SBATCH --output slurm/%a.out

options_file=options/options_run_moments_jackknife.tsv
script_file=run_moments.py

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
