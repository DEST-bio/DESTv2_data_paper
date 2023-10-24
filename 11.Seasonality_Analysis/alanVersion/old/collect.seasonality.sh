#!/usr/bin/env bash
#
#SBATCH -J glmOmn # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00 #<= this may depend on your resources
#SBATCH --mem 15G #<= this may depend on your resources
#SBATCH -p largemem
#SBATCH -A berglandlab_standard
#SBATCH -o /scratch/aob2x/DEST2_analysis/seasonality/logs/collectSeasonality.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST2_analysis/seasonality/logs/collectSeasonality.%A_%a.err # Standard error

###9060

### sbatch --array=1-30 /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/collect.seasonality.sh
###
### sacct -j 50769062
### seff 50160609_376
### cat /scratch/aob2x/DEST2_analysis/seasonality/logs/collectSeasonality.50692434_2.err



module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

cd /home/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion


Rscript \
--vanilla \
collect.seasonality.R \
${SLURM_ARRAY_TASK_ID}
