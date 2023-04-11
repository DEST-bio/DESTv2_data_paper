# ijob -A berglandlab_standard -c20 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])
  model=as.character(args[2])
  #nJobs=as.numeric(args[2])

  ## jobId=1; nJobs=5000

### libraries
  library(data.table)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(5)
  library(tidyverse)

  #setwd("/scratch/aob2x")

### load data
  ### seasonal pairs
    seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
  ### core20  
    core.20 <- fread("./core20_samps.csv")
    names(core.20)[1] = "sampleId_orig"
  ### dest samps  
    samps <- fread("./dest_v2.samps_25Feb2023.csv")
  
    core20.upd = left_join(core.20, samps[,c("sampleId", "sampleId_orig")])

  #### model selector
    if(model == "all_seas") {
      
      message("chosen model --> all")
      seasonal.sets = seasonal.sets
      
    } else if(model == "NoCore20_seas") {
      
      message("chosen model --> No Core20")
      seasonal.sets = seasonal.sets %>% filter(!sampleId %in% core20.upd$sampleId)
      
    } else{ message("model is not specified"); q("no") }
    
    
  ### gds object
    genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds")

  ### sample metadata
  #system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")


### get basic index
  data <- seqGetData(genofile, "annotation/info/AF")
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

### get subset
###
###
  #snp.dt
  snp.dt <- snp.dt[global_af>=0.1]
  #dim(snp.dt)[1] -> total.snps
  #snp.indexes = 1:total.snps
  #colnames(snp.dt) -> SNPguides
  win.bp <- 12000
  step.bp <- 12001
  setkey(snp.dt, "chr")
  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind",
                  .errorhandling="remove")%dopar%{

                    tmp <- snp.dt %>%
                      filter(chr == chr.i)

                    data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }

  wins[,i:=1:dim(wins)[1]]
  ###dim(wins)
  ### ----> 9060

  wins.i = wins[jobId]

  tmp.ids <- snp.dt[chr==wins.i$chr][pos%in%c(wins.i$start:wins.i$end)]$variant.id

  # tmp.ids <- c(1333833, 595225)

### iterate through
  o <- foreach(i=1:length(tmp.ids), .combine="rbind")%do%{
    # i<-1
    message(paste(i, length(tmp.ids), sep=" / "))
    af <- getData(variant=tmp.ids[i])
    af <- merge(af, seasonal.sets, by="sampleId")
    af[,year_pop:=as.factor(interaction(locality, year))]
    af <- af[!is.na(af_nEff)]

    # t1 <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1, data=af, family=binomial())
    # t2 <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season, data=af, family=binomial())
    #
      t3.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ year_pop, data=af,
                     family=quasibinomial())

      t4.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season+year_pop, data=af,
                     family=quasibinomial())

    # t3.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1 + (1|year_pop), data=af,
    #            family=binomial())
    # t4.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season+ (1|year_pop), data=af,
    #            family=binomial())


    obs <-
    data.table(perm=0,
               b_seas=coef(t4.real)[2], se_temp=summary(t4.real)$coef[2,2],
               nTotal=dim(seasonal.sets)[1],
               nObs=dim(af)[1],
               nFixed=sum(af$af_nEff==0) + sum(af$af_nEff==1),
               af=mean(af$af_nEff), neff=mean(af$nEff),
               lrt=anova(t3.real, t4.real, test="Chisq")[2,4],
               p_lrt=anova(t4.real, t3.real, test="Chisq")[2,5],
               model=model
               )



    set.seed(1234)
    nPerm <- 100
    perms <- foreach(j=1:nPerm, .combine="rbind")%dopar%{

      tmp <- af
      tmp[,season:=sample(season)]

      #t3.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~        year_pop,
      #               data=tmp, family=binomial())
      #t4.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season+year_pop,
      #               data=tmp, family=binomial())

      #t3.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1 + (1|year_pop),
      #               data=tmp, family=binomial())
      #t4.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1|year_pop),
      #               data=tmp, family=binomial())

      t3.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ year_pop,
                     data=tmp, family=quasibinomial())
      t4.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + year_pop,
                     data=tmp, family=quasibinomial())


      data.table(perm=j,
                 b_seas=coef(t4.perm)[2], se_temp=summary(t4.perm)$coef[2,2],
                 nTotal=dim(seasonal.sets)[1],
                 nObs=dim(tmp)[1],
                 nFixed=sum(tmp$af_nEff==0) + sum(tmp$af_nEff==1),
                 af=mean(tmp$af_nEff), neff=mean(tmp$nEff),
                 lrt=anova(t3.perm, t4.perm, test="Chisq")[2,4],
                 p_lrt=anova(t4.perm, t3.perm, test="Chisq")[2,5],
                 model=model
      )

    }

    out <- rbind(obs, perms)
    out[,variant.id:=tmp.ids[i]]

  }

  o <- merge(o, snp.dt, by="variant.id")

  #### SAVE O
  output_file = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_ALAN_APR122023/"
  save(o,
       file = paste(output_file,
                    "GLM_out.",
                    jobId,
                    ".",
                    model,
                    ".",
                    paste(wins.i$chr,wins.i$start,wins.i$end, sep = "_"),
                    ".Rdata", sep = ""))
