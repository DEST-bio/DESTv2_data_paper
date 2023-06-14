#!/usr/bin/env bash
#
#SBATCH -J runSNAind # A single job name for the array
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 #<= this may depend on your resources
#SBATCH --mem 120G #<= this may depend on your resources
#SBATCH -o ./slurmOut/snap.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/snap.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH -A jcbnunez


##### More insights into SNAPE
##### This is a bash script
##### JCBN June 8, 2023

module load tabix

## declare SNAPE
snape=/scratch/yey2sn/DEST2_analysis/SNAPE_debugging/snape-pooled/snape-pooled

# Declare parameters

#mpileup="/project/berglandlab/DEST/dest_mapped/RECENT_OUTPUTS/US_Cal_And_-1_2013-10-15/US_Cal_And_-1_2013-10-15.mel_mpileup.txt.gz"
output=./
sample="test_joint"
theta=0.005
D=0.01
priortype="informative"
fold="unfolded"
nflies=100

###

#cp $mpileup ./
#gunzip ./US_Cal_And_-1_2013-10-15.mel_mpileup.txt.gz

awk '{if (last != $1) close(last); print >> $1; last = $1}' US_Cal_And_-1_2013-10-15.mel_mpileup.txt

for chr in {2L,2R,3L,3R,4}; do
  $snape -nchr $(($nflies*2)) -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

echo "done"