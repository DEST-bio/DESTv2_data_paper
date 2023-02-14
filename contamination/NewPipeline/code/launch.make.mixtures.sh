#!/usr/bin/env bash
#
#SBATCH -J sim.read.mix # A single job name for the array
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 #<= this may depend on your resources
#SBATCH --mem 90G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/rmix.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/rmix.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=0-40

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj
#####
Rscript \
--vanilla \
make.sample.mixtures.R \
${SLURM_ARRAY_TASK_ID}

echo "done"