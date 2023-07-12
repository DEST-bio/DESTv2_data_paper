

#module load GCC/10.2.0
#module load Tools/bcftools-1.17
#source activate dest

##created TF with conda with the dependencies
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate tf 

#-----------------------------Launch script----------------------------
# Example
# sbatch --array=1-10 run_relernn2.sh

#-----------------------------Working paths----------------------------
## !!! Modify DIR and DIRDATA: DIR is where ReLERNN2 will write results,
## !!! DIRDATA where the DEST2 vcf is 
DIR="/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/results"
DIRDATA="/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/data"
PBS_ARRAY_INDEX=$(cat /media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/index.txt)

VCFFILE="${DIRDATA}/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz"

#SAMPLE="$(cat $DIR/0-selected_samples.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )"
##SUITABLE FOR OPEN PBS
SAMPLE=$( sed "${PBS_ARRAY_INDEX}q;d" /media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/0-selected_samples.txt )

#-----------------------------Create sample VCF------------------------
VCF="$DIR/VCF"
#mkdir -p $VCF/$SAMPLE
#bcftools view -c1 -Oz -s $SAMPLE -o $VCF/$SAMPLE/$SAMPLE.vcf.gz $VCFFILE
#tabix -p vcf $VCF/$SAMPLE/$SAMPLE.vcf.gz

#-----------------------------Filter VCF--------------------------------
## VCF for autosomes
#bcftools view -Oz -m 2 -M 2 --types=snps --regions 2L,2R,3L,3R $VCF/$SAMPLE/$SAMPLE.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.aut.vcf.gz
#bcftools sort -Oz $VCF/$SAMPLE/$SAMPLE.biallelic_snps.aut.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz
### VCF for X chromsoome
#bcftools view -Oz -m 2 -M 2 --types=snps --regions X $VCF/$SAMPLE/$SAMPLE.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.X.vcf.gz
#bcftools sort -Oz $VCF/$SAMPLE/$SAMPLE.biallelic_snps.X.vcf.gz -o $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz

#-----------------------------Allele freq--------------------------------
AF="$DIR/alleleFreq"
#mkdir -p "$AF/$SAMPLE"
#bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz -i 'AF[0]>0.01 && FORMAT/DP>10' -o $AF/$SAMPLE/$SAMPLE.aut.pool
#bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz -i 'AF[0]>0.01 && FORMAT/DP>10' -o $AF/$SAMPLE/$SAMPLE.X.pool

#-----------------------------Sample depth-------------------------------
DEPTH="$DIR/depth/"
#mkdir -p "$DEPTH/$SAMPLE"
#vcftools --depth --gzvcf $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.aut.vcf.gz --out $DEPTH/$SAMPLE/$SAMPLE.aut
#vcftools --depth --gzvcf $VCF/$SAMPLE/$SAMPLE.biallelic_snps.sort.X.vcf.gz --out $DEPTH/$SAMPLE/$SAMPLE.X

sampleDepthAut=$(cut -f3 $DEPTH/$SAMPLE/$SAMPLE.aut.idepth | tail -n+2 | cut -f1 -d'.')
sampleDepthX=$(cut -f3 $DEPTH/$SAMPLE/$SAMPLE.X.idepth | tail -n+2 | cut -f1 -d'.')

#-----------------------------Number of flies in pool---------------------
nFlies=$(awk  -F,  -v sample="$SAMPLE" ' $1 == sample ' $DIRDATA/dest_v2.samps_8Jun2023.csv | cut -f24 -d',')
nFliesAut=$(( nFlies * 2 ))

#-----------------------------Run ReLERNN---------------------
#mkdir -p ReLERNN_OUTPUT/$SAMPLE/aut
#mkdir -p ReLERNN_OUTPUT/$SAMPLE/X

source ~/miniconda3/etc/profile.d/conda.sh
conda activate relernn2
#conda activate relernn2

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

## Train network
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

## Train network
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
## Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${RELERNNDIRX} \
    --nCPU 32

submitdir="/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2"
cd $submitdir
echo "SUCCESSFUL PROCESS TERMINATION" > controlfile.txt
