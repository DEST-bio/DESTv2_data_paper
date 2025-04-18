#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --time=5:00:00
#SBATCH -o ./SlurmOut/R.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=2-999

module load Rgeospatial

#####
    
   Rscript --vanilla \
     linear_Admix_SLIM.R \
     ${SLURM_ARRAY_TASK_ID}
     
