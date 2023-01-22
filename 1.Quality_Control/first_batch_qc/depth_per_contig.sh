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

# sbatch --array=1-724 ~/DESTv2_data_paper/1.Quality_Control/first_batch_qc/depth_per_contig.sh
# sacct -j 46190508
# cat /scratch/aob2x/dest/slurmOutput/bam_qc.46190508_3.err

### define jobs
  # ls -l /project/berglandlab/DEST/dest_mapped/* | cut -d' ' -f9 > /scratch/aob2x/depth/jobs.csv
  # head -n3 /scratch/aob2x/depth/jobs.csv
  # sed -i '1,3d' /scratch/aob2x/depth/jobs.csv
  # head -n1 /scratch/aob2x/depth/jobs.csv
#
#
  #wc -l /scratch/aob2x/depth/jobs.csv
  #SLURM_ARRAY_TASK_ID=1

  jobN=${SLURM_ARRAY_TASK_ID}
  job=$( sed "${jobN}q;d" /scratch/aob2x/depth/jobs.csv )
  echo $job

  # job=CH_Ors_23_Aug_2020
  ### fo 
### load modules
  module load samtools

### iterate
  cp /project/berglandlab/DEST/*/*/${job}/${job}.original.bam \
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
