#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --time=5:00:00
#SBATCH -o ./SlurmOut/slim.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=2-999

repid=${SLURM_ARRAY_TASK_ID}
root=/users/j/c/jcnunez/scratch/DEST2.0_analysis/REVISION3_ADMIX

mkdir results

module load slim

#####
    
    slim \
	-d "repId=${repid}" \
	-d "root='${root}'" \
     Admix_sim_nonWF.slim
     
mv Agnostic.${SLURM_ARRAY_TASK_ID}.pop.*.16900.txt results
mv AFR.${SLURM_ARRAY_TASK_ID}.pop.*.16900.txt results