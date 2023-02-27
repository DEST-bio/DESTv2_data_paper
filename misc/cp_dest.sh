#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 1G
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch ~/DESTv2_data_paper/misc/cp_dest.sh
### sacct -j 47128513

cp /scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.norep.ann.gds \
/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.5.25Feb2023.norep.ann.gds

cp /scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.norep.ann.vcf.gz \
/project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.5.25Feb2023.norep.ann.vcf.gz

cp /scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.norep.ann.vcf.gz.tbi \
/project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.5.25Feb2023.norep.ann.vcf.gz.tbi
