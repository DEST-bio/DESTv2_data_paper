#!/usr/bin/env bash
#
#SBATCH -J MineGDS # A single job name for the array
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 #<= this may depend on your resources
#SBATCH --mem 120G #<= this may depend on your resources
#SBATCH -p largemem
#SBATCH -A jcbnunez

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

Rscript \
--vanilla \
1.MakePCA.R