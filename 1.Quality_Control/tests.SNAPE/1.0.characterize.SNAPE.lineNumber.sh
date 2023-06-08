#!/usr/bin/env bash
#
#SBATCH -J zcat.runs # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 #<= this may depend on your resources
#SBATCH --mem 9G #<= this may depend on your resources
#SBATCH -o ./slurmOut/zcat.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/zcat.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH --array=1-742

k=${SLURM_ARRAY_TASK_ID}

### Characterize the line numbers of SNAPE outputs in DEST files
### using an array job for efficiency
### June 7, 2023
### Jcbn

# working folder in  RIVANNA --> /scratch/yey2sn/DEST2_analysis/SNAPE_debugging

### Identify the folders to be explored
folders="/project/berglandlab/DEST/dest_mapped"
files=$(find $folders/*/*/*.SNAPE.complete.masked.sync.gz)

### Create file array
myarr=()
for i in $files
do
    myarr+=("$i")
done 

### run zcat command
#for i in ${myarr[1]} ## <== test in one observation
#for 
i=${myarr[k]} ## <== test in five observation
#do
#echo $i
linenum=$(zcat $i | wc -l)
echo -e "$i\t$linenum" >> zcatvalues.DEST.txt

#done


