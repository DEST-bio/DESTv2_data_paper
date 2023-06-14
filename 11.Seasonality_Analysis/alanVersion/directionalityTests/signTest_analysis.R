# ijob -A berglandlab_standard -c5 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)
  library(SeqArray)


### focal
  load("NoCore20_NoProblems_NoFlip_seas_LocRan.Rdata")
  nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
  focal <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]
  rm(mod.out); gc()

### General metadata
  #samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv")


### wd
  setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/glm_output")
  samps[locality=="AT_Wie_Gro"][year==2012]


### seasonal pairs
  seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
  setDT(seasonal.sets)
  seasonal.sets[locality=="AT_Wie_Gro"]

  table(seasonal.sets[,.N,loc.y]$N)
  dim(seasonal.sets)

### add phylo cluster
  phylo_clust <- as.data.table(get(load("~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/phylocluster_data.Rdata")))
  setnames(phylo_clust, "smapleId", "sampleId")

  dim(seasonal.sets)
  seasonal.sets <- merge(seasonal.sets, phylo_clust, by="sampleId", all.x=T)
  dim(seasonal.sets)

  seasonal.sets$cluster[is.na(seasonal.sets$cluster)] = 1

### add in sample metadata
  dim(seasonal.sets)
  seasonal.sets <- merge(seasonal.sets, samps[,c("sampleId", "Recommendation", "exactDate", "continent")], by="sampleId")
  seasonal.sets[is.na(exactDate)]
  dim(seasonal.sets)

  table(seasonal.sets$Recommendation)
  table(seasonal.sets$exactDate)
  table(seasonal.sets$cluster, seasonal.sets$continent)

### trim out ghost samples and leftover pair
  seasonal.sets.ag <- seasonal.sets[,.N,loc.y]
  setkey(seasonal.sets, loc.y)
  dim(seasonal.sets)
  seasonal.sets <- seasonal.sets[J(seasonal.sets.ag[N==2])]
  dim(seasonal.sets)

### target populations
  seasonal.sets <- seasonal.sets

### gds object
  message("open genofile")
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds")

### get basic index
  message("get snp table")
  # data <- seqGetData(genofile, "annotation/info/AF")
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAllele=seqGetData(genofile, "$num_allele"),
                       variant.id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAllele==2]
  seqSetFilter(genofile, snp.dt$variant.id)
  snp.dt[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]

### function
  getData <- function(variant, samps2use=seasonal.sets$sampleId) {
    # variant=snp.dt[chr=="2L"][pos==14617051]$variant.id;
    # samps2use=seasonal.sets$sampleId
    # variant=595225

    ### filter to target
      seqResetFilter(genofile)
      seqSetFilter(genofile, variant.id=variant, sample.id=samps2use)

    ### get frequencies
      message("Allele Freqs")

      ad <- seqGetData(genofile, "annotation/format/AD")$data
      dp <- seqGetData(genofile, "annotation/format/DP")

      af <- data.table(ad=expand.grid(ad)[,1],
                       dp=expand.grid(dp)[,1],
                       sampleId=rep(seqGetData(genofile, "sample.id"),
                                    dim(ad)[2]),
                       variant.id=rep(seqGetData(genofile, "variant.id"),
                                      each=dim(ad)[1]),
                        chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"))

    ### tack them together
      af[,af:=ad/dp]

    ### calculate effective read-depth
      afis <- merge(af, samps[,c("sampleId", "nFlies"), with=F], by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
      afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      afis
  }

### get focal sites

  dat.ag <- foreach(perm.i=c(0:10), .combine="rbind")%dopar%{
    message(perm.i)
    setkey(focal, perm)
    o <- focal[J(perm.i)]
    o[,q:=p.adjust(p_lrt, "fdr")]
    o[,rnp:=rank(p_lrt)/(1+length(p_lrt))]



    rnp.i<-.01
    tmp <- o[rnp<=rnp.i]
    dim(tmp)


    dat <- getData(variant=tmp$variant.id, samps2use=seasonal.sets$sampleId)

    dat <- merge(dat, seasonal.sets, by="sampleId")

    #datw <- dcast(dat, loc.y+variant.id~season, value.var=c("af"))


    dat.ag <- dat[,list(st=af_nEff[season=="spring"] - af_nEff[season=="fall"], nSpring=sum(season=="spring"), nFall=sum(season=="fall"), perm=perm.i, nEff=mean(nEff, na.rm=T), nFlies=mean(nFlies, na.rm=T)),
                   list(loc.y, delta.T.sign, delta.T.mag, T.dir, exactDate, continent, Core20_sat, variant.id)]


    setkey(dat.ag, variant.id, perm)
    setkey(tmp, variant.id, perm)
    merge(dat.ag, tmp)
  }

  save(dat.ag, file="~/signTest.data.Rdata")


### system("scp aob2x@rivanna.hpc.virginia.edu:~/signTest.data.Rdata ~/.")
  library(data.table)
  library(dplyr)
  library(ggplot2)
  load("~/signTest.data.Rdata")
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv")
  samps[,loc.y:=paste(locality, year, sep="_")]
  samps.ag <- samps[,list(lat=mean(lat)), loc.y]
  merge(dat.ag, samps.ag, by="loc.y")

  dat.ag[,focal:=ifelse(delta.T.sign==-1 & delta.T.mag=="Steep" & Core20_sat==F, T, F)]
  dat.ag[,group:=case_when(
        delta.T.sign==-1 & delta.T.mag=="Steep"	& Core20_sat==F  ~ "Focal",
        delta.T.sign==1 & delta.T.mag=="Steep" ~ 	"Steep_Pos",
        Core20_sat==T ~ "Core20",
        TRUE ~ "remainder")]




### SNP-wise analysis
  dat.ag.ag <- dat.ag[,list(pr=mean(st<0, na.rm=T), mu=mean(st, na.rm=T)), list(variant.id, perm, group)]

  ggplot(data=dat.ag.ag, aes(mu)) + geom_histogram() + facet_grid(perm~group)
  ggplot(data=dat.ag.ag, aes(x=pr, y=mu)) + geom_point() + facet_grid(perm~group)

  ggplot(data=dat.ag.ag[perm==0][group!="remainder"], aes(x=group, y=pr, group=variant.id)) + geom_line()

  tmp <- dcast(dat.ag.ag[perm==0][group!="Steep_Pos"], variant.id + perm ~ group, value.var="mu")

  ggplot(data=tmp, aes(x=Focal, y=remainder)) + geom_point()

  fisher.test(table((tmp$Focal)>.0, (tmp$remainder)>.0))

### population wise analysis
  dat.ag.ag <- dat.ag[group=="Focal",list(pr=mean(st<0, na.rm=T), mu=mean(st, na.rm=T)), list(variant.id, perm)]
  setnames(dat.ag.ag, "pr", "focal_pr")

  setkey(dat.ag, variant.id, perm)
  setkey(dat.ag.ag, variant.id, perm)
  dat.ag.pr <- merge(dat.ag, dat.ag.ag)

  bt <- function(x) {

    tab <- table(x)
    message(tab)
    binom.test(tab)$estimate
  }
  dat.ag.pr <- merge(dat.ag.pr, samps.ag, by="loc.y")
    d2 <- dat.ag.pr[,list(TT=sum(sign(st)==sign(mu), na.rm=T),
                          FF=sum(sign(st)!=sign(mu), na.rm=T), nEff=mean(nEff, na.rm=T)),
                    list(loc.y, group, perm, lat)]

    d2[,pr:=TT/(TT+FF)]
    ggplot(data=d2, aes(x=group, y=pr)) + geom_point() + facet_grid(~perm) + ylim(0,1) + geom_hline(yintercept=.5)

    summary(lm(pr~nEff, data=d2[perm==0][group=="Focal"]))


    ggplot(data=d2[group=="Focal"], aes(x=lat, y=pr)) + geom_point() + facet_grid(~perm) + ylim(0,1) + geom_hline(yintercept=.5)


    bto <- foreach(i=1:dim(d2)[1], .combine="rbind", .errorhandling="remove")%do%{
      bt <- binom.test(d2[i]$TT, d2[i]$TT + d2[i]$FF, .5)
      cbind(d2[i], data.table(pval=bt$p.value, est=bt$est, lci=bt$conf.int[1], uci=bt$conf.int[2]))
    }



  ggplot(data=bto, aes(x=group, y=est, color=as.factor(pval<1e-4))) + geom_point() + facet_grid(chr~perm) + ylim(0,1) + geom_hline(yintercept=.5)

  ggplot(data=d2[group=="Focal"], aes(x=T.dir, y=pr)) + geom_point() + facet_grid(~perm)









    dat.ag.pr[perm==0]



    ss.tmp <- seasonal.sets[,list(delta.T.sign=delta.T.sign[1], delta.T.mag=delta.T.mag[1]), loc.y]
    dat.ag <- merge(dat.ag, ss.tmp, by="loc.y")
    dat.ag[,variant.id:=snp.id]
    dat.ag
  }

  m.ag <- m[,list(pr=mean(st<0, na.rm=T)), list(variant.id)]



















  save(o, file="~/NoCore20_NoProblems_NoFlip_seas.perm0.Rdata")

### system("scp aob2x@rivanna.hpc.virginia.edu:~/NoCore20_NoProblems_NoFlip_seas.perm0.Rdata ~/.")

  load("~/NoCore20_NoProblems_NoFlip_seas.perm0.Rdata")

  ggplot(data=o[p_lrt<.01], aes(p_lrt)) + geom_histogram(bins=100)


  gg_qqplot <- function(ps, ci = 0.95) {
    n  <- length(ps)
    df <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
      geom_ribbon(
        mapping = aes(x = expected, ymin = clower, ymax = cupper),
        alpha = 0.1
      ) +
      geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
      # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
      # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
      xlab(log10Pe) +
      ylab(log10Po)
  }

  qqPlot <- gg_qqplot(ps=o$p_lrt)

  ggplot(data=o, aes(x=se_temp, y=(b_seas))) + geom_hex()


  fisher.test(table(o$se_temp<.05, o$p_lrt<.05))

  o[se_temp>.05, q2:=p.adjust(p_lrt, "fdr")]
