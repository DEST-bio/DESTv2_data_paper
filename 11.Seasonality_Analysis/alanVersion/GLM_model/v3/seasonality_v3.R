# ijob -A berglandlab -c5 -p standard --mem=4G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

  args = commandArgs(trailingOnly=TRUE)

  jobId=as.numeric(args[1])
  pops=as.character(args[2]) # all_seas ; NoCore20_seas ; Core20_seas
  nPerm = as.numeric(args[3])
  #model_features=as.character(args[4]) #No_Phylo; Phylo_LocRan; PhyloRan_LocRan; Phylo_Loc; LocRan

  #jobId=450; pops="NoCore20_seas_europe"; nPerm=0

### does file already exist?
  output_dir = "/scratch/aob2x/DEST2_analysis/seasonality/glm_test_May24_2024/"
  output_dir = paste(output_dir, pops, "/", sep="")
  output_dir_final = paste(output_dir, "Loc_Ran_PCA", "/", sep="")
  file = paste(output_dir_final,
               "GLM_out.",
               jobId,
               ".",
               pops,
               ".",
               "omnibus",
               ".",
               "*",
               ".Rdata", sep = "")

  if(length(list.files(output_dir_final, paste("GLM_out.", jobId, ".", pops, ".*.Rdata", sep="")))>0) {
    message("already done")
    q("no")
  }

### libraries
  .libPaths(c("/project/berglandlab/Rlibs_4.3.1/")); .libPaths()

  library(data.table)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(10)
  #library(tidyverse)
  library(lme4)
  library(dplyr)
  # library(ggplot2)
  #setwd("/scratch/aob2x")

### load data

# General metadata
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv")


#### population selector
  if(pops == "all_seas") {

    message("chosen model --> all")
    seasonal.sets = seasonal.sets

  } else if(pops == "NoCore20_seas") {

    message("chosen model --> NoCore20_seas")
    seasonal.sets = seasonal.sets %>% filter(Core20_sat == FALSE)


  } else if(pops == "Core20_seas") {

    message("chosen model --> Only Core20")
    seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE)

  } else if(pops== "NoCore20_NoProblems_seas") {

    message("chosen model --> NoCore20_NoProblems_seas")
    seasonal.sets = seasonal.sets %>% filter(Core20_sat == FALSE) %>% filter(Recommendation == "Pass")
    seasonal.sets.N <- seasonal.sets[,.N,loc.y]
    setkey(seasonal.sets, loc.y)
    seasonal.sets <- seasonal.sets[J(seasonal.sets.N[N==2])]
    dim(seasonal.sets)

  } else if(pops== "NoCore20_NoProblems_Steep_Neg_seas") {

    message("chosen model --> NoCore20_NoProblems_Steep_Neg_seas")
    #seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE) %>% filter(delta.T.sign==-1) %>% filter(delta.T.mag$Steep)
    seasonal.sets <- seasonal.sets[Recommendation == "Pass"][Core20_sat==F][delta.T.sign==-1][delta.T.mag=="Steep"]
    seasonal.sets[,.N,loc.y]

  } else if(pops== "NoCore20_NoProblems_Steep_Pos_seas") {

    message("chosen model --> NoCore20_NoProblems_Steep_Pos_seas")
    #seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE) %>% filter(delta.T.sign==-1) %>% filter(delta.T.mag$Steep)
    seasonal.sets <- seasonal.sets[Recommendation == "Pass"][Core20_sat==F][delta.T.sign==1][delta.T.mag=="Steep"]
    seasonal.sets[,.N,loc.y]

  } else if(pops== "NoProblems_Steep_Pos_seas") {

    message("chosen model --> NoProblems_Steep_Pos_seas")
    #seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE) %>% filter(delta.T.sign==-1) %>% filter(delta.T.mag$Steep)
    seasonal.sets <- seasonal.sets[Recommendation == "Pass"][delta.T.sign==1][delta.T.mag=="Steep"]
    seasonal.sets[,.N,loc.y]

  } else if(pops== "NoCore20_seas_europe") {
    message("chosen model --> NoCore20_seas_europe")
    #seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE) %>% filter(delta.T.sign==-1) %>% filter(delta.T.mag$Steep)

    seasonal.sets <- samps[Recommendation=="Pass"][continent%in%c("Europe", "Asia")][set!="DrosRTEC"][,list(.N,
                                                                   season=c("spring", "fall"),
                                                                   sampleId=c(sampleId[which.min(jday)], sampleId[which.max(jday)]),
                                                                   jday=c(jday[which.min(jday)], jday[which.max(jday)]),
                                                                   diff=c(jday[which.min(jday)]-jday[which.max(jday)])),
                                                                list(country, locality, year, loc_rep, cluster2.0_k5)][N>1][abs(diff)>30]

  } else {
    message("population set is not specified")
    q("no")
  }

### gds object
  message("open genofile")
  #genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds")
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.3May2024.ann.gds")
  #genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.10March2023.ann.gds")


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

      afis[chr=="X", nEff:=round((dp*nFlies)/(dp+nFlies-1))]
      afis[chr!="X", nEff:=round((dp*2*nFlies )/(dp+2*nFlies- 1))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      afis
  }

### make windows and get subset
  message("make windows")
  snp.dt <- snp.dt[global_af>=0.05]
  win.bp <-  50000   #12000
  step.bp <- 50001   #12001
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
  dim(wins)
  ### ----> 2173

  wins.i = wins[jobId]

  tmp.ids <- snp.dt[chr==wins.i$chr][pos%in%c(wins.i$start:wins.i$end)]$variant.id

  if(length(tmp.ids)==0) {
    message("nothing in window")
    q("no")
  }
  # tmp.ids <- tmp.ids[1:3]

### load permutations
  load("/standard/vol186/bergland-lab/alan/dest_baypass/contrast_pop.perm.Rdata")

### iterate through
  message("iterate")
  o <- foreach(i=1:length(tmp.ids), .combine="rbind", .errorhandling="remove")%do%{
    #i <- 1; tmp.ids <- 786996
    message(paste(i, length(tmp.ids), sep=" / "))

    ### get allele frequency data
      af <- getData(variant=tmp.ids[i])   ## af <- getData(variant=2178993)
      af <- merge(af, seasonal.sets, by="sampleId")
      af[,year_pop:=as.factor(interaction(locality, year))]
      af <- af[!is.na(af_nEff)]
      af <- af[af_nEff>0 & af_nEff<1]

    ### iterate through permutations
      o <- foreach(j=0:nPerm, .combine="rbind", .errorhandling="remove")%dopar%{
        # j<-0
          tmp <- af
          cont.perm[perm==j]
          tmp <- merge(tmp, cont.perm[perm==j], by="sampleId")
          setnames(tmp, "season.perm", "season")

        message(j)
        ### iterate through model types

          foreach(model_features = c("yearPop_Ran"),  .combine="rbind", .errorhandling="remove")%do%{
            p_lrt=-999
            seas.AIC = -999
            null.AIC = -999
            # model_features <- "LocBinomial"
            if(model_features == "yearPop_Binomial"){
              message("yearPop_Binomial model")
              t3.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ year_pop,          data = tmp, family= binomial)
              t4.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + year_pop, data = tmp, family= binomial)
              t3.sing <- t4.sing <- NA
            } else if(model_features == "yearPop_Ran" ){
              message("yearPop_Ran model")
              t3.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1      + (1 | year_pop),  data=tmp, family = binomial)
              t4.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1 | year_pop),  data=tmp,  family = binomial)
              t3.sing <- isSingular(t3.real); t4.sing <- isSingular(t4.real)

            }

            if(grepl("Ran", model_features)) {
              p_lrt=anova(t4.real, t3.real, test="Chisq")[2,8]
              seas.AIC = extractAIC(t4.real)[2]
              null.AIC = extractAIC(t3.real)[2]
            } else if (!grepl("Ran", model_features)) {
              p_lrt=anova(t4.real, t3.real, test="Chisq")[2,5]
              seas.AIC = extractAIC(t4.real)[2]
              null.AIC = extractAIC(t3.real)[2]
            }


            obs <-
              data.table(variant.id=tmp.ids[i], perm=j,
                         b_seas=summary(t4.real)$coef[2,1], se_temp=summary(t4.real)$coef[2,2],
                         nTotal=dim(seasonal.sets)[1],
                         nObs=dim(af)[1],
                         nFixed=sum(af$af_nEff==0) + sum(af$af_nEff==1),
                         af=mean(af$af_nEff), neff=mean(af$nEff),
                         p_lrt=p_lrt,
                         pops=pops,
                         model_features=model_features,
                         seas.AIC = seas.AIC,
                         null.AIC = null.AIC,
                         spring.frac = mean(tmp$season=="spring", na.rm=T),
                         ran=runif(1, 0,1e6), singular=paste(t3.sing, t4.sing, sep="_"))
              return(obs)
          } # iterate through models
      } # iterate through perms

      #
      o

    } # last line

  o <- merge(o, snp.dt, by="variant.id")

#### SAVE O
  message("save")
  output_dir = "/scratch/aob2x/DEST2_analysis/seasonality/glm_test_May24_2024/"
  if(!dir.exists(output_dir)) {
    message("makding new dir:")
    message(output_dir)
    dir.create(output_dir)
  }

  output_dir = paste(output_dir, pops, "/", sep="")
  if(!dir.exists(output_dir)) {
    message("makding new dir:")
    message(output_dir)
    dir.create(output_dir)
  }

  foreach(mf=unique(o$model_features))%do%{
    output_dir_final = paste(output_dir, mf, "/", sep="")

    if(!dir.exists(output_dir_final)) {
      message("makding new dir:")
      message(output_dir_final)
      dir.create(output_dir_final)
    }
    oo <- o[model_features==mf]
    save(oo,
         file = paste(output_dir_final,
                      "GLM_out.",
                      jobId,
                      ".",
                      pops,
                      ".",
                      "omnibus",
                      ".",
                      paste(wins.i$chr,wins.i$start,wins.i$end, sep = "_"),
                      ".Rdata", sep = ""))
  }

  message("done")




### old

### seasonal pairs
  # seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
  # seasonal.sets <- get(load("/project/berglandlab/jcbnunez/Shared_w_Alan/seasonal_classification/DEST2.seasonals.plusCore20.flip.met.Rdata"))
  # setDT(seasonal.sets)
  # table(seasonal.sets[,.N,loc.y]$N)
  # dim(seasonal.sets)

### add phylo cluster for 26April2023 version
  #phylo_clust <- as.data.table(get(load("~/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/phylocluster_data.Rdata")))
  #setnames(phylo_clust, "smapleId", "sampleId")
  #
  #dim(seasonal.sets)
  #seasonal.sets <- merge(seasonal.sets, phylo_clust, by="sampleId", all.x=T)
  #dim(seasonal.sets)
  #
  #seasonal.sets$cluster[is.na(seasonal.sets$cluster)] = 1

### add in sample metadata
  # dim(seasonal.sets)
  # seasonal.sets <- merge(seasonal.sets, samps[,c("sampleId", "Recommendation", "exactDate", "continent", "cluster1.0", "cluster2.0_k4",  "cluster2.0_k5", "cluster2.0_k8")], by="sampleId")
  # seasonal.sets[is.na(exactDate)]
  # dim(seasonal.sets)

### set phylo_cluster
  # table(seasonal.sets$Recommendation)
  # table(seasonal.sets$exactDate)
  # table(seasonal.sets$cluster, seasonal.sets$continent)
#
  # table(seasonal.sets$Recommendation, seasonal.sets$cluster2.0_k4)
  # table(seasonal.sets$Recommendation, seasonal.sets$cluster2.0_k5)
  # table(seasonal.sets$Recommendation, seasonal.sets$cluster2.0_k8)

  #seasonal.sets[,cluster:=cluster2.0_k8]
