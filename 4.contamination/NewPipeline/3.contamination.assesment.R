### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(Rsamtools)
library(tidyverse)
library(magrittr)
###
###
####

#folder <- "/gpfs/gpfs0/scratch/yey2sn/DEST2_analysis/kmer.analysis"

####
samps = fread("reads.file.txt", header = F)

### simulans contamination rate
### 
simContam <- foreach(samp.i=samps$V1, .errorhandling="remove")%do%{
  #samp.i=samps$V1[1]
  message(samp.i)
  
  simBam <- paste("map.w.bwamem",samp.i,
                  paste(samp.i,".sim.bam", sep = ""),
                  sep = "/"
                  )
  
  #simBam <- gsub("mark_duplicates_report.txt", "sim.bam", fns[grepl(samp.i, fns)])
  simIdx <- paste(simBam, "bai", sep=".")
  

  melBam <- paste("map.w.bwamem",samp.i,
                  paste(samp.i,".mel.bam", sep = ""),
                  sep = "/"
  )
  
  #simBam <- gsub("mark_duplicates_report.txt", "sim.bam", fns[grepl(samp.i, fns)])
  melIdx <- paste(melBam, "bai", sep=".")
  
  
  
  
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
}

simContam <- rbindlist(simContam)
simContam.ag <- simContam[chr%in%c("2L", "2R", "3L", "3R", "X"), list(propSimNorm=mean(propSimNorm, na.rm=T)), list(auto=chr=="X", sampleId)]

###
simContam.ag %>%
  filter(auto == TRUE) %>%
  separate(sampleId, remove = F, into = c("sim","TrueSim","mel","TrueMel","tot","MixedTot","reads")) %>%
  mutate(TrueSim = as.numeric(TrueSim)/as.numeric(MixedTot) ) %>%
  ggplot(aes(
    x=TrueSim,
    y=propSimNorm,
  )) +
  ylim(0,1) + 
  xlim(0,1) +
  geom_point() +
  geom_abline(slope = 1) ->
  contam.plot

ggsave(contam.plot, file = "contam.plot.pdf")


write.csv(simContam.ag, file="/scratch/aob2x/dest/DEST_freeze1/populationInfo/sequencingStats/simulans.csv", row.names=F)
