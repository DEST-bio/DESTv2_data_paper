#!/usr/bin/env bash
#
#SBATCH -J collect_Core20 # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 02:30:00 #<= this may depend on your resources
#SBATCH --mem 500G #<= this may depend on your resources
#SBATCH -p largemem
#SBATCH -A berglandlab_standard
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.%A_%a.err # Standard error

### sbatch /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/launch.collect.core20.sh
###
### sacct -j 50085870
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/glmOmn.50085870_*.err
### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err



module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1

cd /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion


Rscript \
--vanilla \
collect.seasonality.core20.R
