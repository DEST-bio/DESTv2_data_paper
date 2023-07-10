#!/usr/bin/env bash
#
#SBATCH -J explo # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH -o ./slurmOut/minepca.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/minepca.%A_%a.err # Standard error
#SBATCH --array=1-40 

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

Rscript \
--vanilla \
3.mine.new.pca.loadings.R \
${SLURM_ARRAY_TASK_ID}

