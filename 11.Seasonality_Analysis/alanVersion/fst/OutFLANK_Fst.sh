#!/usr/bin/env bash
#
#SBATCH -J fst # A single job name for the array
#SBATCH -c 20 ### 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:25:00
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/logs/dest_fst.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/dest_fst.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### run as: sbatch --array=1-1086 PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/fst.*.err


### modules
  module load intel/18.0  intelmpi/18.0 R/4.1.1

### run window
  Rscript --vanilla /scratch/aob2x/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/fst/OutFLANK_Fst.R ${SLURM_ARRAY_TASK_ID}
