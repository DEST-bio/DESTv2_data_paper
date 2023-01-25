#!/bin/bash
#
#SBATCH -J runSnakemake # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 01:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/DESTv2_output/logs/runSnakemake.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DESTv2_output/logs/runSnakemake.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

module load samtools/1.9

### sbatch --array=1-687 /scratch/aob2x/DESTv2_data_paper/freebayes/getSmallBams.sh

# find /project/berglandlab/DEST/dest_mapped/ -name '*.mel.bam' > /scratch/aob2x/smallBam/jobs.csv
# SLURM_ARRAY_TASK_ID=1
jobN=${SLURM_ARRAY_TASK_ID}
job=$( sed "${jobN}q;d" /scratch/aob2x/smallBam/jobs.csv )
echo $job

samp=$( echo $job | rev | cut -f1 -d'/' | rev )
echo $samp

region="2L:13613293-13750800"

#samtools index $job

samtools view -b $job $region > /scratch/aob2x/smallBam/small.${samp}
