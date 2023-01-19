#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-62

module load qualimap

JAVAMEM=39G

sample_metadat=../samps.to.map.rd1.txt

sample=$( cat $sample_metadat  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo $sample

file_in=/project/berglandlab/DEST/dest_mapped/pipeline_output/$sample/$sample.mel.bam

qualimap bamqc -bam $file_in \
 -outdir ./$file_in \
 --java-mem-size=$JAVAMEM


date
echo "done"