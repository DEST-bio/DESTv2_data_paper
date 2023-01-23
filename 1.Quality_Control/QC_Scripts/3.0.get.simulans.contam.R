### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(Rsamtools)
library(tidyverse)


### load metadata
samps <- fread("../DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")
fns <- system("ls /project/berglandlab/DEST/dest_mapped/DESTv2/*/*duplicates_report.txt", intern=T)

samps <- filter(samps, set == "DrosEU_3")
#####
### simulans contamination rate
#####

###simContam <- foreach(samp.i=samps[set!="dgn"]$sampleId, .errorhandling="remove")%do%{

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])

samp.i=samps$sampleId[i]

#samp.i=m[set!="dgn"][auto==T]$sampleId[1]
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

idx.out.ag

save(idx.out.ag, file = paste(samp.i, "Rdata", sep ="."))

##}
#simContam <- rbindlist(simContam)
#simContam.ag <- simContam[chr%in%c("2L", "2R", "3L", "3R", "X"), list(propSimNorm=mean(propSimNorm, na.rm=T)), list(auto=chr=="X", sampleId)]
#write.csv(simContam.ag, file="/scratch/aob2x/dest/DEST_freeze1/populationInfo/sequencingStats/simulans.csv", row.names=F)
