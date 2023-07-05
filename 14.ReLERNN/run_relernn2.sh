#!/bin/bash
# -----------------------------Name of the job-------------------------

#SBATCH -J ReLERNN # A single job name for the array
#SBATCH --ntasks=32
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 200G #<= this may depend on your resources
#SBATCH --gres=gpu:p100:1
#SBATCH -p gpu
#SBATCH -A jcbnunez
#SBATCH -o ./slurmOut/ReLERNN.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/ReLERNN.%A_%a.err # Standard error

#-----------------------------Modules----------------------------------
module purge
module load gcc anaconda/2020.11-py3.8

module load bcftools
module load vcftools

source activate relernn2
#-----------------------------Launch script----------------------------
# Example
# sbatch --array=1-10 run_relernn2.sh

#-----------------------------Working paths----------------------------
## !!! Modify DIR and DIRDATA: DIR is where ReLERNN2 will write results,
## !!! DIRDATA where the DEST2 vcf is 
DIR="/scratch/yey2sn/DEST2_analysis/relernn_runs"
DIRDATA="/project/berglandlab/DEST/vcf/"
VCFFILE="${DIRDATA}/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz"
SAMPLE="$(cat $DIR/0-selected_samples.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )"
METADATA="dest_v2.samps_8Jun2023.csv"

echo $SAMPLE
#-----------------------------Create sample VCF------------------------
# This takes a couple minutes

#tabix -p vcf $VCFFILE <-- file must be indexed

VCF="$DIR/VCF"
mkdir -p $VCF/$SAMPLE
bcftools view -c1 -Oz -s $SAMPLE -o $VCF/$SAMPLE/$SAMPLE.vcf.gz $VCFFILE
tabix -p vcf $VCF/$SAMPLE/$SAMPLE.vcf.gz

#-----------------------------Filter VCF--------------------------------
# VCF for autosomes
bcftools view -Oz -m 2 -M 2 --types=snps --regions 2L,2R,3L,3R $VCF/$SAMPLE/$SAMPLE.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.aut.vcf.gz
bcftools sort -Oz $VCF/$SAMPLE/$SAMPLE.biallelic_snps.aut.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz
# VCF for X chromsoome
bcftools view -Oz -m 2 -M 2 --types=snps --regions X $VCF/$SAMPLE/$SAMPLE.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.X.vcf.gz
bcftools sort -Oz $VCF/$SAMPLE/$SAMPLE.biallelic_snps.X.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz

#-----------------------------Allele freq--------------------------------
AF="$DIR/alleleFreq"
mkdir -p "$AF/$SAMPLE"
bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz -i 'AF[0]>0.01 && FORMAT/DP>10' -o $AF/$SAMPLE/$SAMPLE.aut.pool
bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz -i 'AF[0]>0.01 && FORMAT/DP>10' -o $AF/$SAMPLE/$SAMPLE.X.pool

#-----------------------------Sample depth-------------------------------
DEPTH="$DIR/depth/"
mkdir -p "$DEPTH/$SAMPLE"
vcftools --depth --gzvcf $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz --out $DEPTH/$SAMPLE/$SAMPLE.aut
vcftools --depth --gzvcf $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz --out $DEPTH/$SAMPLE/$SAMPLE.X

sampleDepthAut=$(cut -f3 $DEPTH/$SAMPLE/$SAMPLE.aut.idepth | tail -n+2 | cut -f1 -d'.')
sampleDepthX=$(cut -f3 $DEPTH/$SAMPLE/$SAMPLE.X.idepth | tail -n+2 | cut -f1 -d'.')

echo $sampleDepthAut
echo $sampleDepthX
#-----------------------------Number of flies in pool---------------------
nFlies=$(awk  -F,  -v sample="$SAMPLE" ' $1 == sample ' $METADATA | cut -f24 -d',')
nFliesAut=$(( nFlies * 2 ))

echo $nFlies
echo $nFliesAut

#-----------------------------Run ReLERNN---------------------
mkdir -p ReLERNN_OUTPUT/$SAMPLE/aut
mkdir -p ReLERNN_OUTPUT/$SAMPLE/X

#conda activate relernn2 #<-- this happens above

# ReLERNN parameters
SIMULATE="ReLERNN_SIMULATE_POOL"
TRAIN="ReLERNN_TRAIN_POOL"
PREDICT="ReLERNN_PREDICT_POOL"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="1"
MU="3.27e-9"    # Assumed per-base mutation rate
GENTIME="0.08"  # Assumed generation time (in years)
URTR="10"
RELERNNDIRAUT="$DIR/ReLERNN_OUTPUT/$SAMPLE/aut"
RELERNNDIRX="$DIR/ReLERNN_OUTPUT/$SAMPLE/X"
POOL_AUT="$AF/$SAMPLE/$SAMPLE.aut.pool"
POOL_X="$AF/$SAMPLE/$SAMPLE.X.pool"
GENOME_AUT="$DIR/genome.bed"
GENOME_X="$DIR/X.bed"

# FOR AUTOSOMES
# Simulate data
${SIMULATE} \
    --pool ${POOL_AUT} \
    --sampleDepth ${nFliesAut} \
    --genome ${GENOME_AUT} \
    --projectDir ${RELERNNDIRAUT} \
    --assumedMu ${MU} \
    --assumedGenTime ${GENTIME} \
    --upperRhoThetaRatio ${URTR} \
    --nCPU 32

# Train network
${TRAIN} \
    --projectDir ${RELERNNDIRAUT} \
    --readDepth ${sampleDepthAut} \
    --maf 0.01 \
    --nCPU 32

# Predict
${PREDICT} \
    --pool ${POOL_AUT} \
    --projectDir ${RELERNNDIRAUT} \
    --minSites 50

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${RELERNNDIRAUT} \
    --nCPU 32

# FOR THE X
# Simulate data
${SIMULATE} \
    --pool ${POOL_X} \
    --sampleDepth ${nFlies} \
    --genome ${GENOME_X} \
    --projectDir ${RELERNNDIRX} \
    --assumedMu ${MU} \
    --assumedGenTime ${GENTIME} \
    --upperRhoThetaRatio ${URTR} \
    --nCPU 32

# Train network
${TRAIN} \
    --projectDir ${RELERNNDIRX} \
    --readDepth ${sampleDepthX} \
    --maf 0.01 \
    --nCPU 32

# Predict
${PREDICT} \
    --pool ${POOL_X} \
    --projectDir ${RELERNNDIRX} \
    --minSites 50

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${RELERNNDIRX} \
    --nCPU 32
    
rm -rf ${RELERNNDIRAUT}/train

rm -rf ${RELERNNDIRX}/train
rm -rf $VCF/$SAMPLE
