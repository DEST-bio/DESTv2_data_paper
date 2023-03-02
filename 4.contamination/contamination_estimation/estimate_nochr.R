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
                list(nSim=as.numeric(tstrsplit(samp, "\\.")[[2]]), nMel=as.numeric(tstrsplit(samp, "\\.")[[4]]))]

  sim.idx.ag[,simNorm:=sim_mapped/chrLen_sim]
  sim.idx.ag[,melNorm:=mel_mapped/chrLen_mel]

  sim.idx.ag[,simRate_norm:=simNorm/(simNorm+melNorm)]

  sim.idx.ag[,simRate:=sim_mapped/(sim_mapped+mel_mapped)]


  #t1 <- lm(I(nSim/(nSim+nMel))~simRate_norm*chr + I(simRate_norm^2)*chr, sim.idx.ag)
  t1 <- lm(I(nSim/(nSim+nMel))~simRate + I(simRate^2), sim.idx.ag)

  sim.pred <- as.data.table(expand.grid(simRate=seq(from=0, to=1, by=.01)))
  sim.pred[,pred:=predict(t1, sim.pred)]

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
                list(samp=tstrsplit(samp, "/")[[3]])]


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

  #real.idx.ag[!is.na(chr),auto:=chr%in%c("2L", "2R", "3L", "3R")]

save(real.idx.ag, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/contamination_estimation/simulans_rates.Rdata")
### filter

  model_plot <-
  ggplot(data=sim.pred, aes(x=simRate, y=pred)) +
  geom_abline(slope=1, intercept=0) +
  geom_line() +
  geom_point(data=sim.idx.ag, aes(x=simRate, y=I(nSim/(nSim+nMel)))) +
  xlab("Observed simulans contamination\nbased on mapping rate to hologenome\n(simRate)") + ylab("\'True\' simulans contamination rate\n(simulated)")


  C <- ggplot(data=real.idx.ag, aes(x=simRate, y=simRate_kmer)) + geom_point() + geom_abline(intercept=0, slope=1) + ylim(0,1) + xlim(0,1)

  model_plot + C


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

  distributionA <- ggplot(data=real.idx.ag, aes(pred_simRate)) + geom_histogram(bins=20) +
                  facet_grid(~set) +
                  scale_y_continuous(trans = "mylog10") +
                  theme(axis.text.x=element_text(size=8, angle=90))

  distributionB <- ggplot(data=real.idx.ag, aes(simRate_kmer)) + geom_histogram(bins=20) +
                  facet_grid(~set) +
                  scale_y_continuous(trans = "mylog10") +
                  theme(axis.text.x=element_text(size=8, angle=90))


  layout <- "
  AB
  CC
  DD"

  mega <- model_plot + C + distributionA + distributionB + plot_layout(design=layout) + plot_annotation(tag_level="A")
  ggsave(mega, file="~/contamination.pdf")



  layout <- "
  AAABCD
  AAABCD
  EEEBCD
  EEEBCD"

  model_plot + A + B + C + distribution + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')
