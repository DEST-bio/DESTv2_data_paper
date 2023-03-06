#!/bin/bash
#
#SBATCH -J redovcf # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 ### 1 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/DESTv2_output/logs/redovcf.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DESTv2_output/logs/redovcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch ~/DESTv2_data_paper/misc/remakevcf/run_redovcf.sh
### sacct -j 47379537

module load parallel

popSet=all; method=PoolSNP; maf="001"; mac=50; version=25Feb2023; wd=/scratch/aob2x/DESTv2_outputPoolSNP50


cd /home/aob2x/DESTv2_data_paper/misc/remakevcf

# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 2L
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 2R
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 3L
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 3R
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 X
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 Y
# ./gather.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 4
#
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 2L
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 2R
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 3L
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 3R
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 X
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 Y
# ./gather.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 4

# ./annotate.sh all PoolSNP "001" 50 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 ~/snpEff
./annotate.sh PoolSeq SNAPE "NA" NA 25Feb2023 /scratch/aob2x/DESTv2_outputPoolSNP50 ~/snpEff
