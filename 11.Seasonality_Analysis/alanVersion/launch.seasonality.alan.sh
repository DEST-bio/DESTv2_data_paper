#!/usr/bin/env bash
#
#SBATCH -J glmAlLo # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A berglandlab_standard
#SBATCH -o ./slurmOut/glmAlLo.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/glmAlLo.%A_%a.err # Standard error
#SBATCH --array=1-1000

###9060

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

Rscript \
--vanilla \
seasonality.R \
${SLURM_ARRAY_TASK_ID} \
"all_seas" \
"Phylo_LocRan"


#Rscript \
#--vanilla \
#seasonality.R \
#${SLURM_ARRAY_TASK_ID} \
#"NoCore20_seas"