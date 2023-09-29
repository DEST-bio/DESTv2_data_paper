#!/usr/bin/env bash
#
#SBATCH -J runSNAind # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 36:00:00 #<= this may depend on your resources
#SBATCH --mem 60G #<= this may depend on your resources
#SBATCH -o ./slurmOut/snap.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/snap.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH -A jcbnunez
#SBATCH --array=1-16


##### More insights into SNAPE
##### This is a bash script
##### JCBN June 8, 2023

module load tabix

## declare SNAPE -- software
snape=/scratch/yey2sn/DEST2_analysis/Old_Scratch/SNAPE_debugging/snape-pooled/snape-pooled


#### declare samples
sample_file=/scratch/yey2sn/DEST2_analysis/New_Scratch/individual_snape/collapsed.samples.txt

rootfile=/project/berglandlab/DEST/dest_mapped/RECENT_OUTPUTS

### extract samples
sampleid=$( cat $sample_file  | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo $sampleid

# Declare parameters

mpileup=$rootfile/$sampleid/$sampleid.mel_mpileup.txt.gz
echo $mpileup

####
mkdir $sampleid
cd $sampleid
cp $mpileup ./

## unzip
gunzip $sampleid.mel_mpileup.txt.gz

#####
output=./
sample=$sampleid
theta=0.005
D=0.01
priortype="informative"
fold="unfolded"
nflies=100

###
awk '{if (last != $1) close(last); print >> $1; last = $1}' $sampleid.mel_mpileup.txt

for chr in {2L,2R,3L,3R,4}; do
  $snape -nchr $(($nflies*2)) -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

echo "done"