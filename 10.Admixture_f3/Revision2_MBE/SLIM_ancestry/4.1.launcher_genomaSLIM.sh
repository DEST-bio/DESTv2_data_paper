#!/bin/sh
#
#SBATCH -J moments_admix
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-999

module load Rgeospatial

Rscript --vanilla 4.genomalicous.SliM.R \
${SLURM_ARRAY_TASK_ID} 
