#### ---> Re run the contamnation assesment for DEST 2.0
#### Collapsed samples

### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(Rsamtools)
library(tidyverse)
library(vroom)

### --> -bash-4.2$pwd /scratch/yey2sn/DEST2_analysis/filtering
#outfile = "/scratch/yey2sn/DEST2_analysis/filtering/sim_contam_final/"
outfile = "/gpfs/gpfs0/scratch/yey2sn/DEST2_analysis/filtering/collapsed_samps"  

### load metadata
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.qc_merge.csv"

samps <- fread(meta_git)
setDT(samps)
#samps = samps[set!="dgn"]
samps = samps[!is.na(collapsedSamples)]

fns <- system("ls /project/berglandlab/DEST/dest_mapped/*/*/*duplicates_report.txt", intern=T)
#fns <- system("ls /project/berglandlab/DEST/dest_mapped/DEST_plus_collapsed/*/*duplicates_report.txt", intern=T)

##log.dat = vroom("./log.simcontam.txt",  col_names = F)
##names(log.dat) = c("sample.ith", "status")
##process.samples =
##foreach(i=1:711, .combine = "rbind")%do%{
##  data.frame(sample.ith = str_split(fns, "/")[[i]][7])
##}
##left_join(process.samples, log.dat) %>%
##  mutate(i = 1:dim(.)[1]) %>%
##  filter(is.na(status)) %>%
##  .$i

###
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
samp.i=samps$sampleId[i]
#

#simContam <- foreach(samp.i=samps[set!="dgn"]$sampleId, .errorhandling="remove")%do%{
#simContam <- foreach(samp.i=missingsamps, .errorhandling="remove")%do%{

#samp.i=samps[set!="dgn"]$sampleId[1]

message(samp.i)
simBam <- gsub("mark_duplicates_report.txt", "sim.bam", fns[grepl(samp.i, fns)])
simIdx <- paste(simBam, "bai", sep=".")

melBam <- gsub("mark_duplicates_report.txt", "mel.bam", fns[grepl(samp.i, fns)])
melIdx <- paste(melBam, "bai", sep=".")

system( paste("module load samtools; samtools index ", simBam) )
system( paste("module load samtools; samtools index ", melBam) )

simidx.out <- as.data.table(idxstatsBam(file=simBam, index=simIdx))[grepl("2L|2R|3L|3R|X|^4$|Y", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]
melidx.out <- as.data.table(idxstatsBam(file=melBam, index=melIdx))[grepl("2L|2R|3L|3R|X|^4$|Y", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]

idx.out <- merge(melidx.out, simidx.out, by="seqnames")

idx.out[,nReads:=mapped.x + mapped.y]
idx.out[,simChr:=grepl("sim", seqnames)]
idx.out[,chr:=gsub("sim_", "", seqnames)]
idx.out[,nReadsNorm:=nReads/seqlength.x]

idx.out.ag <- idx.out[,list(propSim=nReads[simChr==T]/sum(nReads),
                            propSimNorm=nReadsNorm[simChr==T]/sum(nReadsNorm),
                            nMelReads=nReads[simChr==F],
                            melChrLen=seqlength.x[simChr==F],
                            sampleId=samp.i),
                      list(chr)]
idx.out.ag[,mappingEffort:=nMelReads/melChrLen]

#idx.out.ag

save(idx.out.ag, 
     file = paste(outfile,"SimContam.",samp.i, ".Rdata", sep ="")
)

#} ## close loop
###

##  write.table(data.frame(s=samp.i, i="Complete"), 
##              file = "log.simcontam.txt", 
##              append = TRUE, quote = FALSE, sep = "\t",
##              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
##              col.names = FALSE, 
##              qmethod = c("escape", "double"),
##              fileEncoding = "")

