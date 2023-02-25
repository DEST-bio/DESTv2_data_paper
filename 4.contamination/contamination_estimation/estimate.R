### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)
  library(readxl)

### simulated data to build model
  load("~/simContam.Rdata")

  sim.dat <- simContam
  sim.dat[,chr:=tstrsplit(V1, "\t")[[1]]]
  sim.dat[,chrLen:=as.numeric(tstrsplit(V1, "\t")[[2]])]
  sim.dat[,nReads:=as.numeric(tstrsplit(V1, "\t")[[3]])]

  sim.dat[grepl("sim_", chr), species:="simulans"]
  sim.dat[grepl("Scaffold", chr), species:="melanogaster"]
  setkey(sim.dat, chr)
  sim.dat[J(c("2L", "2R", "3L", "3R", "4", "X", "Y")), species:="melanogaster"]

  sim.dat[J(c(paste("sim_", c("2L", "2R", "3L", "3R", "X"), sep=""), c("2L", "2R", "3L", "3R" ,"X"))), muller:=T]
  sim.dat[is.na(muller), muller:=F]
  sim.dat[J(c(paste("sim_", c("2L", "2R", "3L", "3R", "X"), sep=""), c("2L", "2R", "3L", "3R" ,"X"))), chrSimp:=gsub("sim_", "", chr)]

  sim.idx.ag <-
           sim.dat[,list(nReads=sum(nReads, na.rm=T),
                     chrLen_mel=sum(chrLen[species=="melanogaster"], na.rm=T),
                     chrLen_sim=sum(chrLen[species=="simulans"], na.rm=T),
                     mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                     muller="muller"),
                list(chr=chrSimp, nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))]

  sim.idx.ag[,simNorm:=sim_mapped/chrLen_sim]
  sim.idx.ag[,melNorm:=mel_mapped/chrLen_mel]

  sim.idx.ag[,simRate_norm:=simNorm/(simNorm+melNorm)]

  sim.idx.ag[,simRate:=sim_mapped/(sim_mapped+mel_mapped)]


  t1 <- lm(I(nSim/(nSim+nMel))~simRate_norm*chr + I(simRate_norm^2)*chr, sim.idx.ag)

  sim.pred <- as.data.table(expand.grid(simRate_norm=seq(from=0, to=1, by=.01), chr=c("2L", "2R", "3L", "3R", "X")))
  sim.pred[,pred:=predict(t1, sim.pred)]

  model_plot <-
  ggplot(data=sim.pred, aes(x=simRate_norm, y=pred, group=chr, color=chr)) +
  geom_abline(slope=1, intercept=0) +
  geom_line() +
  geom_point(data=sim.idx.ag, aes(x=simRate_norm, y=I(nSim/(nSim+nMel)), group=chr, color=chr)) +
  xlab("Observed simulans contamination based on reads") + ylab("\'True\' simulans contamination rate")

### real data
  load("~/simContam_destv2.Rdata")

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

  real.idx.ag <-
           dat[,list(nReads=sum(nReads, na.rm=T),
                     chrLen_mel=sum(chrLen[species=="melanogaster"], na.rm=T),
                     chrLen_sim=sum(chrLen[species=="simulans"], na.rm=T),
                     mel_mapped=sum(nReads[species=="melanogaster"], na.rm=T),
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T),
                     muller="muller"),
                list(chr=chrSimp, samp=tstrsplit(samp, "/")[[3]])]


  real.idx.ag[,samp:=gsub(".original.bam", "", samp)]

  real.idx.ag[,simNorm:=sim_mapped/chrLen_sim]
  real.idx.ag[,melNorm:=mel_mapped/chrLen_mel]

  real.idx.ag[,simRate_norm:=simNorm/(simNorm+melNorm)]
  real.idx.ag[,simRate:=sim_mapped/(sim_mapped+mel_mapped)]

### add in kmer
  kmer <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/readDepth_vs_kmer_v3/dest_kmer_screening.xls"))
  setnames(kmer, "File", "samp")
  setnames(kmer, "Dsimu/(Dmela+Dsimu)", "simRate_kmer")

  real.idx.ag <- merge(real.idx.ag, kmer, all=T)
  table(is.na(real.idx.ag$simRate))
  table(is.na(real.idx.ag$simRate_kmer))

### add in sample metadata
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_21Feb2023.csv")

  real.idx.ag <- merge(real.idx.ag, samps, by.x="samp", by.y="sampleId", all=T)

### predict contamination rate based on model
  real.idx.ag[,pred_simRate:=predict(t1, real.idx.ag)]

  real.idx.ag[,nSim:=nFlies*simRate]

  real.idx.ag[!is.na(chr),auto:=chr%in%c("2L", "2R", "3L", "3R")]



### filter
  ggplot(data=real.idx.ag, aes(x=chr, y=log10(mel_mapped/chrLen_mel), group=samp)) + geom_line() + facet_grid(~set)




  real.idx.ag.ag <- real.idx.ag[,
                            list(nSim=mean(nSim,na.rm=T), simRate_norm=mean(simRate_norm, na.rm=T), normReadDepth=mean(log10(mel_mapped/chrLen_mel)), nFlies=mean(nFlies)), list(auto, samp, set)]



  B <- ggplot(data=real.idx.ag.ag[nSim<3][set!="DEST_plus"][set!="dgn"][!is.na(auto)], aes(y=nSim, x=auto, group=samp)) +geom_line() + facet_grid(set~.) + ggtitle("~")

  C <- ggplot(data=real.idx.ag.ag[nSim<3][set!="DEST_plus"][set!="dgn"][auto==T][!is.na(auto)], aes(nSim)) +
  geom_vline(xintercept=c(.5, 1)*2, color="red") +
  geom_histogram(bins=100) +
  coord_flip() + xlim(-.2, 3) + facet_grid(set~.) + ggtitle("Autosome")

  A <- ggplot(data=real.idx.ag.ag[nSim<3][set!="DEST_plus"][set!="dgn"][auto==F][!is.na(auto)], aes(nSim)) +geom_histogram(bins=100) +
  geom_vline(xintercept=c(.5, 1)*2, color="red") +
  coord_flip() + scale_y_reverse()+ facet_grid(set~.) + ggtitle("X")


  library(scales)

  mylog10_trans <- function (base = 10)
  {
    trans <- function(x) log(x + 1, base)
    inv <- function(x) base^x
    trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
              domain = c(1e-100, Inf))
  }

  ggplot(df, aes(x=x)) +
    geom_histogram() +
    scale_y_continuous(trans = "mylog10")

  distribution <- ggplot(data=real.idx.ag.ag[auto==T], aes(simRate_norm)) + geom_histogram(bins=100) + facet_wrap(~set) +   scale_y_continuous(trans = "mylog10")




  layout <- "
  AAABCD
  AAABCD
  EEEBCD
  EEEBCD"

  model_plot + A + B + C + distribution + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')
