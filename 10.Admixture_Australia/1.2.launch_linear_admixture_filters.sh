#!/usr/bin/env bash
#
#SBATCH -J linearadx # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH -o ./slurmOut/admx.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/admx.%A_%a.err # Standard error
#SBATCH --array=1-347

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

## All was already ran...
for filter in chr2L chr2R chr3L chr3R noINV 2Lt 2RNS 3RK 3RK_3RP_3RMo silent
do
echo $filter

#Now run the script
Rscript \
--vanilla \
1.2.per.CHR.perINV.linear_admixture.R \
${SLURM_ARRAY_TASK_ID} \
$filter

done



