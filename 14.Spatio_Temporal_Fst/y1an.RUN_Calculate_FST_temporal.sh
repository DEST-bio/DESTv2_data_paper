#!/usr/bin/env bash
#
#SBATCH -J fst.time.win # A single job name for the array
#SBATCH -c 30
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 40G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fstim.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/fstim.%A_%a.err # Standard error
#SBATCH --array=1-923

module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

###
n="${SLURM_ARRAY_TASK_ID}"

Rscript \
--vanilla \
y1an.Calculate_array_FST_temporal.oneyear.R ${n}

date
echo "done"