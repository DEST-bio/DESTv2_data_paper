#!/bin/sh
#
#SBATCH -J moments_admix
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-114

module load Rgeospatial

guide=/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/pairs_guide_file.txt

Rscript --vanilla 1.Genomalicious.R \
${SLURM_ARRAY_TASK_ID} \
$guide
