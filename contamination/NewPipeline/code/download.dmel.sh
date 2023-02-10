#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 7:00:00 ### 7 hours
#SBATCH --mem 40G
#SBATCH -o ../slurmOutput/sra.fetch.%A_%a.out # Standard output
#SBATCH -e ../slurmOutput/sra.fetch.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-80

module load sratoolkit/2.10.5

sranum=$(sed '1d' Dmel.80.txt | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F ' ' '{print $2}' )

echo  $sranum

fastq-dump \
--outdir ./ \
${sranum}

gzip ${sranum}.fastq

