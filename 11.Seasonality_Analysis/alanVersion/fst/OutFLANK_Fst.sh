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

### run as: sbatch --array=1-215 ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/fst/OutFLANK_Fst.sh
### sacct -j 54723734
### cat /scratch/aob2x/logs/dest_fst.54723734_136.err


### modules
  module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1
  #module load intel/18.0  intelmpi/18.0 R/4.1.1

### run window
  #SLURM_ARRAY_TASK_ID=1
  Rscript --vanilla ~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/fst/OutFLANK_Fst.R ${SLURM_ARRAY_TASK_ID}
