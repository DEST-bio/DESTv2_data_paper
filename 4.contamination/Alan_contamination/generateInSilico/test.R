# ijob -A berglandlab -c4 -p dev --mem=10G
# module load gcc/9.2.0  openmpi/3.1.6 R/4.2.1 samtools/1.12; R

### libraries
  library(data.table)
  library(foreach)

###
  setwd("/project/berglandlab/DEST/Contamination_test/mapping_output")
  files <- system("ls -d */*.srt.flt.bam", intern=T)

  simContam <- foreach(samp.i=files, .errorhandling="remove")%do%{
    #samp.i=files[1]
    message(samp.i)

    dat <- as.data.table(system(paste("samtools idxstats ", samp.i, sep=""), intern=T))
    dat[,samp:=samp.i]
    dat
  }
  simContam <- rbindlist(simContam)

  save(simContam, file="~/simContam.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/simContam.Rdata ~/.


### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)

### data
  load("~/simContam.Rdata")

  dat <- simContam
  dat[,chr:=tstrsplit(V1, "\t")[[1]]]
  dat[,chrLen:=as.numeric(tstrsplit(V1, "\t")[[2]])]
  dat[,nReads:=as.numeric(tstrsplit(V1, "\t")[[3]])]

  dat[grepl("sim_", chr), species:="simulans"]
  dat[grepl("Scaffold", chr), species:="melanogaster"]
  setkey(dat, chr)
  dat[J(c("2L", "2R", "3L", "3R", "4", "X", "Y")), species:="melanogaster"]

  dat[J(c(paste("sim_", c("2L", "2R", "3L", "3R", "X"), sep=""), c("2L", "2R", "3L", "3R" ,"X"))), muller:=T]
  dat[is.na(muller), muller:=F]

  idx.ag <- rbind(
            dat[muller==T,list(nReads=sum(nReads, na.rm=T),
                      mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                      sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                      muller="muller"),
                 list(nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))],
           dat[,list(nReads=sum(nReads, na.rm=T),
                     mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                     muller="all"),
                list(nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))]
            )

  idx.ag[,simRate:=sim_mapped/nReads]

   mel <- ggplot(data=idx.ag, aes(x=nSim, y=sim_mapped, color=muller)) + geom_line() + ggtitle("Reads mapping to Dsim") + ylim(0,50000)
   sim <- ggplot(data=idx.ag, aes(x=nSim, y=mel_mapped, color=muller)) + geom_line() + ggtitle("Reads mapping to Dmel") + ylim(0,50000)
   tot <- ggplot(data=idx.ag, aes(x=nSim, y=nReads,     color=muller)) + geom_line() + ggtitle("Reads mapping") +             ylim(0,50000)

   idx.ag[order(nSim)]

   corPlot <- ggplot(data=idx.ag, aes(x=nSim/(nSim+nMel), y=simRate, color=muller)) + geom_line() + geom_abline()

   layout <- "
   AB
   CD"

   mel + sim + tot + corPlot +plot_layout(design=layout)



   summary(lm(simRate~I(nSim/(nSim+nMel)), idx.ag[muller=="all"]))
