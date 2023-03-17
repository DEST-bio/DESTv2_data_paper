### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)
  library(readxl)

### simulated data to build model - read depth method
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
                     sim_mapped=sum(nReads[species=="simulans"], na.rm=T)),
                list(nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))]

  sim.idx.ag[,simNorm:=sim_mapped/chrLen_sim]
  sim.idx.ag[,melNorm:=mel_mapped/chrLen_mel]

  sim.idx.ag[,simRate_norm:=simNorm/(simNorm+melNorm)]

  sim.idx.ag[,simRate:=sim_mapped/(sim_mapped+mel_mapped)]

  sim.idx.ag[,simRate2:=simRate^2]
  t1 <- lm(I(nSim/(nSim+nMel))~simRate , sim.idx.ag)
  t1a <- lm(I(nSim/(nSim+nMel))~simRate_norm + I(simRate_norm^2), sim.idx.ag)

  sim.pred <- as.data.table(expand.grid(simRate=seq(from=0, to=1, by=.01)))
  sim.pred[,simRate2:=simRate^2]
  sim.pred[,pred:=predict(t1, sim.pred)]

### load in kmer results
  kmersim <- fread("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/contamination_estimation/tmp_est.csv")
  setnames(kmersim, "Nmin5_conf0.95", "kmerEst")
  kmersim[,kmerEst:=kmerEst/100]
  kmersim[,ExpectedPercentage:=ExpectedPercentage/100]

  t2 <- lm(ExpectedPercentage~kmerEst , kmersim)

### combine
  sim_comb <- rbind(data.table(simRate=sim.idx.ag$simRate, nSim=sim.idx.ag$nSim, method="readDepth"),
                    data.table(simRate=kmersim$kmerEst, nSim=kmersim$NindDsimulans, method="kmer"))




#kmer_model_plot <- ggplot(data=kmersim, aes(x=kmerEst, y=ExpectedPercentage)) +
#geom_abline(slope=1, intercept=0) +
#geom_point(color="red") +
#geom_abline(slope=coef(t2)[2], intercept=coef(t2)[1], color="red")

#model_plot <-
#ggplot(data=sim.pred, aes(x=simRate, y=pred)) +
#geom_abline(slope=1, intercept=0) +
#geom_line() +
#geom_point(data=sim.idx.ag, aes(x=simRate_norm, y=I(nSim/(nSim+nMel)))) +
#xlab("Observed simulans contamination based on reads") + ylab("\'True\' simulans contamination rate")

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
  kmer <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/contamination_estimation/dest_kmer_screening_v2.xls"))
  setnames(kmer, "File", "samp")
  setnames(kmer, "100*Dsimu/(Dmela+Dsimu)", "simRate_kmer")

  real.idx.ag <- merge(real.idx.ag, kmer, all=T)
  table(is.na(real.idx.ag$simRate))
  table(is.na(real.idx.ag$simRate_kmer))

### add in sample metadata
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_21Feb2023.csv")

  real.idx.ag <- merge(real.idx.ag, samps, by.x="samp", by.y="sampleId", all=T)

### predict contamination rate based on model
  real.idx.ag[,pred_simRate:=(coef(t1)[1] + coef(t1)[2]*simRate)]
  real.idx.ag[,pred_kmer:=(coef(t2)[1] + coef(t2)[2]*simRate_kmer)]

  real.idx.ag[,nSim:=nFlies*simRate]
  real.idx.ag[,nSim_kmer:=nFlies*pred_kmer/100]

  real.idx.ag[!is.na(chr),auto:=chr%in%c("2L", "2R", "3L", "3R")]


### filter
  real.idx.ag.ag <- real.idx.ag[,
                            list(nSim=mean(nSim,na.rm=T), nSim_kmer=mean(nSim_kmer,na.rm=T),
                                  simRate=mean(simRate, na.rm=T),
                                  simRate_kmer=mean(simRate_kmer, na.rm=T)/100,
                                  normReadDepth=mean(log10(mel_mapped/chrLen_mel)),
                                  nFlies=mean(nFlies)),
                            list(auto, samp, set)]


  real_contam_plot <- ggplot(data=real.idx.ag.ag[auto==T], aes(x=simRate, y=simRate_kmer)) + geom_point() +
  xlab("simulans rate, read depth") +
  ylab("simunans rate, kmer")


  popCounts <- real.idx.ag.ag[auto==T, list(.N), list(set, simRate_bin=round(simRate, 2))]

  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

  hist_plot <- ggplot(data=popCounts) +
  geom_segment(aes(x=simRate_bin, xend=simRate_bin, yend=0, y=N)) +
  facet_wrap(~set, nrow=1) +
  scale_y_continuous(trans = "log10", breaks = breaks, minor_breaks = minor_breaks) +
  annotation_logticks(outside=T, sides="l") +
  ylab("simulans rate, read_depth")

  nSimPlot <- ggplot(data=real.idx.ag.ag[nSim<3][auto==T], aes(nSim)) + geom_histogram(bins=100) + xlab("Predicted number of contaminants")


   model_plot <- ggplot(data=sim_comb, aes(x=simRate, y=nSim/80, group=method, color=method)) +
   geom_abline(slope=1, intercept=0) +
   geom_point() + geom_line() + ylab("True Dsim proportion") + xlab("Estimated Dsim Proportion") +
   theme(
        legend.position=c(.2,.8)
        )



  layout <- "
  ABC
  DDD"

  mega <- model_plot + real_contam_plot  + nSimPlot + hist_plot +
  plot_layout(design=layout) +
  plot_annotation(tag_levels="A")

ggsave(mega, file="~/contamination.pdf", h=8, w=12)






    B <- ggplot(data=real.idx.ag.ag[nSim<3][set!="DEST_plus"][set!="dgn"][!is.na(auto)], aes(y=nSim_kmer/100, x=auto, group=samp)) +geom_line() + facet_grid(set~.) + ggtitle("~")

    C <- ggplot(data=real.idx.ag.ag[nSim<3][set!="DEST_plus"][set!="dgn"][auto==T][!is.na(auto)], aes(nSim_kmer/100)) +
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

    distribution <- ggplot(data=real.idx.ag.ag[auto==T], aes(simRate)) + geom_histogram(bins=100) + facet_wrap(~set) +
    scale_y_continuous(trans = "mylog10")




  layout <- "
  AABBCDE
  AABBCDE
  FFFFCDE
  FFFFCDE"

  model_plot + kmer_model_plot + A + B + C + distribution + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')
