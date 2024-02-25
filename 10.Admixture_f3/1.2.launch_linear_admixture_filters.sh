#!/usr/bin/env bash
#
#SBATCH -J linearadx # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 40G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fstim.%A_%a.out # Standard output
#SBATCH --array=1-347

module load Rtidyverse

## All was already ran...
# silent
for filter in chr2L chr2R chr3L chr3R noINV 2Lt 2RNS 3RK 3RK_3RP_3RMo All
do
echo $filter

#Now run the script
Rscript \
--vanilla \
1.2.per.CHR.perINV.linear_admixture.R \
${SLURM_ARRAY_TASK_ID} \
$filter

done



