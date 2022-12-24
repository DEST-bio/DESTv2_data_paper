#!/usr/bin/env bash
#
#SBATCH -J bam_qc # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ###
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/bam_qc.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

#ijob -c20 --mem=20G -p standard -A berglandlab

# sbatch --array=1-61 ~/home/aob2x/DESTv2_data_paper/first_batch_qc/depth_per_contig.sh


### define jobs
  # ls -l /project/berglandlab/DEST/dest_mapped/DESTv2/ | cut -d' ' -f9 > /scratch/aob2x/depth/jobs.csv
  # ls -l /project/berglandlab/DEST/dest_mapped/pipeline_output | cut -d' ' -f9 >> /scratch/aob2x/depth/jobs.csv
  # sed -i '/^$/d' /scratch/aob2x/depth/jobs.csv

  jobN=${1}
  job=$( sed "${jobN}q;d" /scratch/aob2x/depth/jobs.csv )
  echo $job

### load modules
  module load samtools

### iterate
  cp /project/berglandlab/DEST/dest_mapped/DESTv2/${job}/${job}.original.bam \
  /scratch/aob2x/depth/.

  samtools index -@ 20 /scratch/aob2x/depth/${job}.original.bam

  ~/mosdepth.1 \
  -t 20 \
  -n \
  -x \
  /scratch/aob2x/depth/${job} \
  /scratch/aob2x/depth/${job}.original.bam

  rm /scratch/aob2x/depth/${job}.original.bam
  rm /scratch/aob2x/depth/${job}.original.bam.bai
