#!/usr/bin/env bash
#
#SBATCH -J nasa # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 00:10:00 ###
#SBATCH --mem 2G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab_standard

# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


# sbatch --array=1-745%10 ~/DESTv2_data_paper/seasonal/getNASApower.sh
# sacct -j 45972442
# cat /scratch/aob2x/dest/slurmOutput/bam_qc.45972442_1.err

# SLURM_ARRAY_TASK_ID=1

Rscript ~/DESTv2_data_paper/seasonal/getNASApower.R ${SLURM_ARRAY_TASK_ID}
