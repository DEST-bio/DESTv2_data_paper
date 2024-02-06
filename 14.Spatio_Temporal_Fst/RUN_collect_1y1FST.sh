#!/usr/bin/env bash
#
#SBATCH -J collweah # A single job name for the array
#SBATCH -c 30
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 #<= this may depend on your resources
#SBATCH --mem 40G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/collweah.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/collweah.%A_%a.err # Standard error

module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

###
Rscript \
--vanilla \
collect.outputs.1ymode.R

date
echo "done"