#!/usr/bin/env bash
#
#SBATCH -J LO.ma # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/win.%A_%a.out # Standard output
#SBATCH --array=2
##3-931

module load Rtidyverse

k="${SLURM_ARRAY_TASK_ID}"

echo $n

Rscript \
--vanilla \
Complete_lift_over_Machado.R $k

date
echo "done"