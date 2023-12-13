### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(scales)

### simulated data to build model
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/data/simContam.Rdata") ### from remapping read counts of simulated data; script?

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

### kmer simulations
  n.ind.simu=seq(0,80,2)
  n.ind.mela=80-n.ind.simu

  nkmin=c(1,5) ; conf=c(0.90,0.95)
  cdt=expand.grid(nkmin,conf)

  res=array(0,dim=c(length(n.ind.simu),length(nkmin)*length(conf),6),
              dimnames = list(n.ind.simu,paste0("Nmin",cdt[,1]," conf",cdt[,2]),
            c("NseqSimu","NseqMela","NseqAnan","NseqWolb","NseqOthers","NseqAssigned")))
  setwd("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/data/ana_simu")
  for(i in 1:length(n.ind.simu)){
    for(j in 1:nrow(cdt)){
      tmp.f=paste0("res_files/sim.",n.ind.simu[i],"_nkmin",cdt[j,1],"_conf",cdt[j,2],".summary")
      tmp.d=read.table(tmp.f,skip=1,row.names=1)
      tmp.v=c(tmp.d["Dsimu",1],tmp.d["Dmela",1],tmp.d["Danan",1],tmp.d["Wolb",1])
      res[i,j,]=c(tmp.v,sum(tmp.d[,1])-sum(tmp.v),sum(tmp.d[,1]))
    }
  }
  tmp.d=100*res[,,1]/(res[,,1]+res[,,2])
  rmse=sqrt(colSums((tmp.d/100-n.ind.simu/80)**2)/nrow(tmp.d))
  kmer.sim <- as.data.table(tmp.d)
  kmer.sim[,nSim_individuals:=n.ind.simu]
  setnames(kmer.sim, names(kmer.sim), gsub(" ", "_", names(kmer.sim)))
  setnames(kmer.sim, names(kmer.sim), gsub("\\.", "", names(kmer.sim)))
  kmer.sim[,simRate:=nSim_individuals/80]
  kmer.sim[,method:="kmer"]
  kmer.sim[,simRate_norm:=Nmin5_conf095/100]

  kmer.sim[,c("simRate_norm", "method", "simRate"), with=F]
  sim.idx.ag[,method:="samtools"]
  sim.idx.ag <- rbind(kmer.sim, sim.idx.ag, fill=T)

### real data

  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/data/simContam_destv2.Rdata") ### made by '/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/old/readDepth_vs_kmer_v3/readDepth_contamination.R'

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


### get out total read counts
  dat.totals <- dat[,list(n_total=sum(nReads, na.rm=T)), list(samp)]
  dat.totals[,samp:=tstrsplit(samp, "/")[[3]]]
  dat.totals[,samp:=gsub(".original.bam", "", samp)]

  real.idx.ag <- merge(real.idx.ag, dat.totals, by="samp")

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
  real.idx.ag[,pred_simRate:=predict(t1, real.idx.ag)]
  real.idx.ag[,nSim:=nFlies*simRate]

### save
  save(real.idx.ag, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/4.contamination/contamination_estimation/simulans_rates.Rdata")


### kmer long
  kl <- melt(kmer[,-c("Dmela+Dsimu+Wolb", "simRate_kmer"), with=F], id.vars="samp")
  kl.ag <- kl[,list(max_rate=max(value), mean_rate=mean(value)), list(variable)]
  kl.ag[,r:=rank(-max_rate)]
  kl.ag <- kl.ag[order(-r)]
  kl.ag[,sp.factor:=factor(variable, levels=kl.ag$variable)]
  kl <- merge(kl, kl.ag, by="variable")
  kl <- merge(kl, samps[,c("sampleId", "nFlies", "set"), with=F], by.x="samp", by.y="sampleId")


### plot
  mylog10_trans <- function (base = 10) {
    trans <- function(x) log(x + 1, base)
    inv <- function(x) base^x
    trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base),
              domain = c(1e-100, Inf))
  }

  model_plot <-
  ggplot(data=sim.idx.ag, aes(x=simRate, y=simRate_norm, group=method, color=method)) +
  geom_abline(slope=1, intercept=0) +
  geom_line() +
  geom_point() +
  xlab("Estimated simulans contamination") + ylab("\'True\' simulans contamination rate\n(simulated)")

  C <- ggplot(data=real.idx.ag, aes(x=simRate, y=simRate_kmer/100)) + geom_point() + geom_abline(intercept=0, slope=1) + ylim(0,1) + xlim(0,1)

  distributionA <- ggplot(data=real.idx.ag, aes(pred_simRate)) + geom_histogram(bins=20) +
                  facet_grid(~set) +
                  scale_y_continuous(trans = "mylog10") +
                  theme(axis.text.x=element_text(size=8, angle=90))

  distributionB <- ggplot(data=real.idx.ag, aes(simRate_kmer)) + geom_histogram(bins=20) +
                  facet_grid(~set) +
                  scale_y_continuous(trans = "mylog10") +
                  theme(axis.text.x=element_text(size=8, angle=90))

  wacky <-   ggplot(data=kl[max_rate>1], aes(x=sp.factor, y=plogis(value/100 * nFlies), group=samp)) +
    geom_point() + geom_line() + facet_grid(~set) + coord_flip()

  layout <- "
  AB
  AB
  CC
  DD"

  mega <- model_plot + C + distributionA + wacky + plot_layout(design=layout) + plot_annotation(tag_level="A")
  ggsave(mega, file="~/contamination.pdf", height=7, width=8)



  layout <- "
  AAABCD
  AAABCD
  EEEBCD
  EEEBCD"

  model_plot + A + B + C + distribution + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')


#### stas
  sim.idx.ag[,list(rmse=sqrt(mean(abs(simRate_norm - simRate), na.rm=T)), cor=cor(simRate_norm, simRate)), list(method)]




  model_plot + C

  ggplot(df, aes(x=x)) +
    geom_histogram() +
    scale_y_continuous(trans = "mylog10")







#### old

### combine with wolbachia estimates
  load("~/wolb.Rdata")
  wolb <- merge(wolb, dat.totals, by.x="sampleId", by.y="samp")
  wolb[,FPK:=nReads/(chrLen/1e3)]
  wolb[,FPKM:=FPK/(n_total/1e6)]
  wolb.wide <- dcast(wolb, sampleId~chr, value.var="FPKM")
  wolb.l <- merge(wolb, real.idx.ag, by.x="sampleId", by.y="samp")

  ggplot(data=wolb.l, aes(x=simRate_kmer, y=FPKM, color=chr)) + geom_point() + facet_grid(chr~set)
  ggplot(data=wolb.l, aes(x=Wolb, y=FPKM, color=chr)) + geom_point() + facet_grid(chr~set)
  ggplot(data=wolb.wide, aes(x=W_pipientis, y=sim_wMel)) + geom_point()

  kmer.samp <- merge(kmer, samps[,c("sampleId", "nFlies", "set"), with=F], by.x="samp", by.y="sampleId")

  ggplot(data=kmer.samp, aes(x=Dmela, y=Wolb)) + geom_point() + facet_grid(~set)


  scale_x_discrete(labels=kl.ag[max_rate>1][order(max_rate)]$variable) +

  ggplot(data=kl.ag, aes(x=max_rate, y=mean_rate)) + geom_point()
