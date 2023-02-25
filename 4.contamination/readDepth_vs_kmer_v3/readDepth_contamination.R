# ijob -A berglandlab -c20 -p standard --mem=40G
# module load gcc/9.2.0  openmpi/3.1.6 R/4.2.1 samtools/1.12; R

### libraries
  library(data.table)
  library(foreach)
  library(foreach)
  library(doMC)
  registerDoMC(20)

###
  setwd("/project/berglandlab/DEST/dest_mapped")
  files <- system("ls -d */*/*.original.bam", intern=T)
  files <- files[!grepl("RECENT_OUTPUTS", files)]

  simContam <- foreach(samp.i=files, .errorhandling="remove")%dopar%{
    #samp.i=files[1]
    message(samp.i)

    if(length(system(paste("ls -d ", samp.i, ".bai", sep=""), intern=T))==0) {
      system(paste("samtools index ", samp.i, sep=""))
    }

    dat <- as.data.table(system(paste("samtools idxstats ", samp.i, sep=""), intern=T))
    dat[,samp:=samp.i]
    dat
  }
  simContam <- rbindlist(simContam)

  save(simContam, file="~/simContam_destv2.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/simContam_destv2.Rdata ~/.


### libraries
  library(ggpubr)
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)
  library(readxl)

### data
  load("~/simContam_destv2.Rdata")

  dat <- simContam
  dat[,chr:=tstrsplit(V1, "\t")[[1]]]
  dat[,chrLen:=as.numeric(tstrsplit(V1, "\t")[[2]])]
  dat[,nReads:=as.numeric(tstrsplit(V1, "\t")[[3]])]

  dat[grepl("sim_", chr), species:="simulans"]
  dat[grepl("Scaffold", chr), species:="melanogaster"]
  setkey(dat, chr)
  dat[J(c("2L", "2R", "3L", "3R", "4", "X", "Y")), species:="melanogaster"]



  real.idx.ag <- rbind(
            dat[muller==T,list(nReads=sum(nReads, na.rm=T),
                      mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                      sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                      muller="muller"),
                 list(samp=tstrsplit(samp, "/")[[3]])],
           dat[,list(nReads=sum(nReads, na.rm=T),
                     mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                     muller="all"),
                list(samp=tstrsplit(samp, "/")[[3]])],
          dat[X==T,list(nReads=sum(nReads, na.rm=T),
                    mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                    sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                    muller="X"),
               list(samp=tstrsplit(samp, "/")[[3]])],
         dat[Y==T,list(nReads=sum(nReads, na.rm=T),
                   mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                   sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                   muller="Y"),
              list(samp=tstrsplit(samp, "/")[[3]])]

            )

  real.idx.ag[,simRate:=sim_mapped/nReads]
  real.idx.ag[,samp:=gsub(".original.bam", "", samp)]

  ggplot(data=real.idx.ag, aes(simRate)) + geom_histogram() + facet_grid(~muller)

  table(real.idx.ag$simRate>.5, real.idx.ag$muller)


### load in kmer
  kmer <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/readDepth_vs_kmer_v3/dest_kmer_screening.xls"))
  setnames(kmer, "File", "samp")
  setnames(kmer, "Dsimu/(Dmela+Dsimu)", "simRate_kmer")

  m <- merge(real.idx.ag, kmer, all=T)
  table(is.na(m$simRate))
  table(is.na(m$simRate_kmer))

### samples
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_21Feb2023.csv")

  n <- merge(m, samps, by.x="samp", by.y="sampleId", all=T)

### load in simulation
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
  dat[J(c(paste("sim_", c("2L", "2R", "3L", "3R", "X"), sep=""), c("2L", "2R", "3L", "3R" ,"X"))), chrSimp:=gsub("sim_", "", chr)]

  sim.idx.ag <- rbind(
            dat[muller==T,list(nReads=sum(nReads, na.rm=T),
                      chrLen_mel=sum(chrLen[species=="melanogaster"], na.rm=T),
                      chrLen_sim=sum(chrLen[species=="simulans"], na.rm=T),
                      mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                      sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                      muller="muller"),
                 list(chrSimp, nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))],
           dat[,list(nReads=sum(nReads, na.rm=T),
                     chrLen_mel=sum(chrLen[species=="melanogaster"], na.rm=T),
                     chrLen_sim=sum(chrLen[species=="simulans"], na.rm=T),
                     mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                     muller="all"),
                list(chrSimp, nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))]
            )
  sim.idx.ag[,simNorm:=sim_mapped/chrLen_sim]
  sim.idx.ag[,melNorm:=mel_mapped/chrLen_mel]

  sim.idx.ag[,simRate:=simNorm/(simNorm+melNorm)]

  t1 <- lm(I(nSim/(nSim+nMel))~simRate*chrSimp + I(simRate^2)*chrSimp, sim.idx.ag[muller=="all"])


  ggplot(data=sim.idx.ag[muller=="all"], aes(x=I(nSim/(nSim+nMel)), y=simRate, group=chrSimp, color=chrSimp)) + geom_point() +
  geom_abline(intercept=0, slope=1) + ylim(0,1) + xlim(0,1) +
  stat_cor(method = "pearson", label.x = .05, label.y = .95) + ggtitle("simulated simulas rate vs. read depth simrate")




  n[,pred_simRate:=predict(t1, n[,"simRate", with=F])]
  n[,nSim:=pred_simRate*nFlies]


  n.ag <- n[,list(diff=nSim[muller=="X"] - nSim[muller=="all"]), list(samp)]
  hist(n.ag$diff, breaks=100)



  A <-
  ggplot(data=n[muller=="all"], aes(x=simRate, y=simRate_kmer)) + geom_point() +
  geom_abline(intercept=0, slope=1) + ylim(0,1) + xlim(0,1) +
  stat_cor(method = "pearson", label.x = .05, label.y = .95) + ggtitle("kmer simuland rate vs. read depth simulans rate")

  B <-
  ggplot(data=idx.ag[muller=="all"], aes(x=I(nSim/(nSim+nMel)), y=simRate)) + geom_point() +
  geom_abline(intercept=0, slope=1) + ylim(0,1) + xlim(0,1) +
  stat_cor(method = "pearson", label.x = .05, label.y = .95) + ggtitle("simulated simulas rate vs. read depth simrate")

  C <-
  ggplot(data=n[muller=="all"][nSim<10], aes(nSim)) + geom_histogram(bins=100) +
  geom_vline(xintercept=c(.5), color="red") +
  geom_vline(xintercept=c(1), color="blue")  +
  ggtitle("Predicted number of simulans")


  ggplot(data=n[muller%in%c("all", "X")][nSim<5], aes(x=muller, y=nSim, group=samp)) + geom_line() +
  geom_hline(yintercept=c(.5), color="red") +
  geom_hline(yintercept=c(1), color="blue")



  A + B + C



### read depth
  dat[J(c(paste("sim_", c("2L", "2R", "3L", "3R", "X"), sep=""), c("2L", "2R", "3L", "3R" ,"X"))), chrSimp:=gsub("sim_", "", chr)]
  dat.ag <- dat[,list(rd_norm=nReads/chrLen), list(species, chr=chrSimp, samp=tstrsplit(samp, "/")[[3]])]
  dat.ag <- na.omit(dat.ag)

  dat.ag.ag <- dat.ag[,list(XA=rd_norm[chr=="X"]/median(rd_norm[chr!="X"])), list(species, samp)]





  dat.ag.ag[,samp:=gsub(".original.bam", "", samp)]
  #dat.ag.ag <- merge(dat.ag.ag, n, by.x="samp", by.y="sampleId")
  dat.ag.ag <- merge(dat.ag.ag, n[muller=="all"], by="samp")


  ggplot(data=dat.ag.ag, aes(XA, group=species, fill=species)) + geom_histogram(bins=100) + facet_grid(set~species) + geom_vline(xintercept=.5)
  ggplot(data=dat.ag.ag[nSim<5], aes(y=XA, x=nSim)) + geom_point() + facet_grid(set~species) + geom_vline(xintercept=.5)


  ggplot(dat.ag, aes(x=chr, y=rd_norm, group=samp)) + geom_line() + facet_grid(~species)






  dat.ag <- dat[!is.na(chr),list(rd_norm_mel=mean(nReads[species=="melanogaster"]/chrLen[species=="melanogaster"]),
                                 rd_norm_sim=mean(nReads[species=="simulans"]/chrLen[species=="simulans"])),
                 list(chr=chrSimp, samp=tstrsplit(samp, "/")[[3]])]


   dat.ag <- dat[!is.na(chr),list(rd_norm_mel=mean(nReads/chrLen),
                                  rd_norm_sim=mean(nReads/chrLen)),
                  list(species, chr=chrSimp, samp=tstrsplit(samp, "/")[[3]])]



  dat.ag <- na.omit(dat.ag)


  ggplot(data=dat.ag, aes(x=rd_norm_mel, y=species, group=samp)) + facet_grid(~chr) + geom_line()
