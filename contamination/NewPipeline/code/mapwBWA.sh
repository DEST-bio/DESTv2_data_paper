#!/usr/bin/env bash
#
#SBATCH -J mapmix # A single job name for the array
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 #<= this may depend on your resources
#SBATCH --mem 90G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/mapmix.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/mapmix.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-41


#Load Rivanna modules 
module load gcc/9.2.0
module load bwa/0.7.17
module load bbmap
module load fastqc
module load samtools
module load qualimap
module load picard

###############

#cp /project/berglandlab/DEST/paramTest/holo_dmel_6.12.tar.gz ./
tar -xvf holo_dmel_6.12.tar.gz

####
samplei=$( cat reads.file.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )


##
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory
REFERENCE=./dmel_holo_sim.fa
###
WORKING_FOLDER=/scratch/yey2sn/DEST2_analysis/kmer.analysis/map.w.bwamem
###

############

mkdir $WORKING_FOLDER/$samplei

###########################################################################
	
bwa mem \
-M \
-t $CPU \
$REFERENCE \
$samplei.fq \
> $WORKING_FOLDER/$samplei/$samplei.sam
		
###########################################################################

# Sort merge bams
java -Xmx$JAVAMEM \
 -jar $PICARD SortSam \
 I=$WORKING_FOLDER/$samplei/$samplei.sam \
 O=$WORKING_FOLDER/$samplei/$samplei.srt.bam \
 SO=coordinate \
 VALIDATION_STRINGENCY=SILENT

# Remove duplicates of final file
java -Xmx$JAVAMEM \
 -jar $PICARD MarkDuplicates \
 I=$WORKING_FOLDER/$samplei/$samplei.srt.bam  \
 O=$WORKING_FOLDER/$samplei/$samplei.srt.rmdp.bam  \
 M=$WORKING_FOLDER/$samplei/$samplei.dupstat.txt \
 VALIDATION_STRINGENCY=SILENT \
 REMOVE_DUPLICATES=true

###
#Filter out the simulans contaminants
mel_chromosomes="2L 2R 3L 3R 4 X Y mitochondrion_genome"
sim_chromosomes="sim_2L sim_2R sim_3L sim_3R sim_4 sim_X sim_mtDNA"

#create species specific bamfiles
samtools view -@ $CPU $WORKING_FOLDER/$samplei/$samplei.srt.rmdp.bam $mel_chromosomes -b > $WORKING_FOLDER/$samplei/$samplei.mel.bam
samtools view -@ $CPU $WORKING_FOLDER/$samplei/$samplei.srt.rmdp.bam $sim_chromosomes -b > $WORKING_FOLDER/$samplei/$samplei.sim.bam


 