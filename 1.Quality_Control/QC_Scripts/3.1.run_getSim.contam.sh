#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 30:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 40G
#SBATCH -o ./slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e ./slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-253

module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

### SLURM_ARRAY_TASK_ID=1
  Rscript ./get.simulans.contam.R ${SLURM_ARRAY_TASK_ID} 

