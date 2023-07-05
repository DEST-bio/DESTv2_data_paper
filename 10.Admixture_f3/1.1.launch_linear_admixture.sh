#!/usr/bin/env bash
#
#SBATCH -J linearadx # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH -o ./slurmOut/admx.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/admx.%A_%a.err # Standard error
#SBATCH --array=1-347

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

Rscript \
--vanilla \
1.0.linear_admixture.R ${SLURM_ARRAY_TASK_ID}