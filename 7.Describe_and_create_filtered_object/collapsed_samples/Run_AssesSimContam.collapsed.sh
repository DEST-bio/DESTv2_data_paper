#!/usr/bin/env bash
#
#SBATCH -J simcontam # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 30:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 90G
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-16

module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

### SLURM_ARRAY_TASK_ID=1
  Rscript ./contam_collapsedsamps.R ${SLURM_ARRAY_TASK_ID} 

