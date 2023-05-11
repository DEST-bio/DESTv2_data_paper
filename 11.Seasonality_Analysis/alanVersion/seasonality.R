# ijob -A berglandlab_standard -c20 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

  args = commandArgs(trailingOnly=TRUE)
  
  jobId=as.numeric(args[1])
  model=as.character(args[2]) # all_seas ; NoCore20_seas ; Core20_seas
  nPerm = as.numeric(args[3])
  #model_features=as.character(args[4]) #No_Phylo; Phylo_LocRan; PhyloRan_LocRan; Phylo_Loc; LocRan
  
  #jobId=1
  #model="all_seas"
  #model_features="Phylo_LocRan"
    
### libraries
  library(data.table)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(5)
  library(tidyverse)
  library(lme4)
  
  #setwd("/scratch/aob2x")

### load data

  # General metadata
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv")  
  
  ### seasonal pairs
  seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
    setDT(seasonal.sets)
    
    phylo_clust <- get(load("phylocluster_data.Rdata"))
    
    left_join(seasonal.sets, phylo_clust ) ->
      seasonal.sets
    
    #seasonal.sets %>% filter(is.na(cluster))
    seasonal.sets$cluster[which(is.na(seasonal.sets$cluster))] = 1
    
  ### core20  
  #  core.20 <- fread("./core20_samps.csv")
  #  names(core.20)[1] = "sampleId_orig"
  #### dest samps  
  #  samps <- fread("./dest_v2.samps_25Feb2023.csv")
  #
  #  core20.upd = left_join(core.20, samps[,c("sampleId", "sampleId_orig")])
  #  seasonal.phylo.clusters = get(load("phylocluster_data.Rdata"))
    
  #### model selector
    if(model == "all_seas") {
      
      message("chosen model --> all")
      seasonal.sets = seasonal.sets
      
    } else if(model == "NoCore20_seas") {
      
      message("chosen model --> No Core20")
      seasonal.sets = seasonal.sets %>% filter(Core20_sat == FALSE)
      
      
    } else if(model == "Core20_seas") {
      
      message("chosen model --> Only Core20")
      seasonal.sets = seasonal.sets %>% filter(Core20_sat == TRUE)
      
    } else{ message("model is not specified"); q("no") }
    
  ### gds object
    genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds")

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

  # tmp.ids <- c(759953, 1333833, 595225)

### iterate through
o.mods = foreach(model_features = c("No_Phylo", 
                                    "Phylo_LocRan", 
                                    "PhyloRan_LocRan", 
                                    "Phylo_Loc", 
                                    "LocRan" ##,  
                                    ), 
                 .combine = "rbind")%do%{

  o <- foreach(i=1:length(tmp.ids), .combine="rbind")%do%{
    
    message(paste(i, length(tmp.ids), sep=" / "))
    
    af <- getData(variant=tmp.ids[i])
    af <- merge(af, seasonal.sets, by="sampleId")
    af[,year_pop:=as.factor(interaction(locality, year))]
    af <- af[!is.na(af_nEff)]
    
    left_join(af, phylo_clust) -> af
    af$cluster = as.factor(af$cluster)
    
    #No_Phylo; Phylo_LocRan
    if(model_features == "No_Phylo"){
      
      message("No_Phylo model")
      t3.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ year_pop,
                     data=af, family=binomial)
      t4.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + year_pop,
                     data=af, family=binomial)
      
      #lrt=anova(t3.real, t4.real, test="Chisq")$`Chisq`[2],
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chi)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    } else if(model_features == "Phylo_LocRan" ){
      
      message("Phylo_LocRan model")
      t3.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster + (1 | year_pop),
                       data = af, family = binomial)
      t4.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster + (1 | year_pop),
                       data = af, family = binomial)
      
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chisq)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    } else if(model_features == "PhyloRan_LocRan" ){
      
      message("PhyloRan_LocRan model")
      t3.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ (1 | cluster) + (1 | year_pop),
                       data = af, family = binomial)
      t4.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1 | cluster) + (1 | year_pop),
                       data = af, family = binomial)
      
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chisq)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    } else if(model_features == "Phylo_Loc" ){
      
      message("Phylo_Loc model")
      t3.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster + year_pop,
                       data = af, family = binomial)
      t4.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster + year_pop,
                       data = af, family = binomial)
      
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chi)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    } else if(model_features == "LocRan" ){
      
      message("LocRan model")
      t3.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ (1 | year_pop),
                       data = af, family = binomial)
      t4.real <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1 | year_pop),
                       data = af, family = binomial)
      
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chisq)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    } else if(model_features == "JustPhylo" ){
      
      message("JustPhylo model")
      t3.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster,
                     data = af, family = binomial)
      t4.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster,
                     data = af, family = binomial)
      
      p_lrt=anova(t4.real, t3.real, test="Chisq")$`Pr(>Chi)`[2]
      seas.AIC = extractAIC(t4.real)[2]
      null.AIC = extractAIC(t3.real)[2]
      
    }
    
    
    
    obs <-
      data.table(perm=0,
                 b_seas=summary(t4.real)$coef[2,1], se_temp=summary(t4.real)$coef[2,2],
                 nTotal=dim(seasonal.sets)[1],
                 nObs=dim(af)[1],
                 nFixed=sum(af$af_nEff==0) + sum(af$af_nEff==1),
                 af=mean(af$af_nEff), neff=mean(af$nEff),
                 #lrt=anova(t3.real, t4.real, test="Chisq")$`Chisq`[2],
                 p_lrt=p_lrt,
                 model=model,
                 model_features=model_features,
                 seas.AIC = seas.AIC,
                 null.AIC = null.AIC
      )
    
    
    set.seed(1234)
    nPerm <- nPerm
    perms <- foreach(j=1:nPerm, .combine="rbind")%dopar%{
      
      tmp <- af
      tmp[,season:=sample(season)]
  
      #No_Phylo; Phylo_LocRan
      if(model_features == "No_Phylo" ){
        
        t3.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ year_pop,
                       data=tmp, family=binomial)
        t4.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + year_pop,
                       data=tmp, family=binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chi)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
      } else if(model_features == "Phylo_LocRan" ){
        
        t3.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster + (1 | year_pop),
                         data = tmp, family = binomial)
        t4.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster + (1 | year_pop),
                         data = tmp, family = binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chisq)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
        
      } else if(model_features == "PhyloRan_LocRan" ){
        
        t3.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ (1 | cluster) + (1 | year_pop),
                         data = tmp, family = binomial)
        t4.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1 | cluster) + (1 | year_pop),
                         data = tmp, family = binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chisq)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
      } else if(model_features == "Phylo_Loc" ){
        
        t3.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster + year_pop,
                         data = tmp, family = binomial)
        t4.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster + year_pop,
                         data = tmp, family = binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chi)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
      }  else if(model_features == "LocRan" ){
        
        t3.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ (1 | year_pop),
                         data = tmp, family = binomial)
        t4.perm <- glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + (1 | year_pop),
                         data = tmp, family = binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chisq)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
      } else if(model_features == "JustPhylo" ){
        
        t3.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ cluster,
                       data = tmp, family = binomial)
        t4.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ season + cluster,
                       data = tmp, family = binomial)
        
        p_lrt=anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chi)`[2]
        seas.AIC = extractAIC(t4.perm)[2]
        null.AIC = extractAIC(t3.perm)[2]
        
      }
      
      
      
      data.table(perm=j,
                 b_seas=summary(t4.perm)$coef[2,1], se_temp=summary(t4.perm)$coef[2,2],
                 nTotal=dim(seasonal.sets)[1],
                 nObs=dim(tmp)[1],
                 nFixed=sum(tmp$af_nEff==0) + sum(tmp$af_nEff==1),
                 af=mean(tmp$af_nEff), neff=mean(tmp$nEff),
                 #lrt=anova(t3.perm, t4.perm, test="Chisq")$`Chisq`[2],
                 p_lrt= p_lrt, #anova(t4.perm, t3.perm, test="Chisq")$`Pr(>Chisq)`[2],
                 model=model,
                 model_features=model_features,
                 seas.AIC = seas.AIC,
                 null.AIC = null.AIC
      )
      
    }
    
    out <- rbind(obs, perms)
    out[,variant.id:=tmp.ids[i]]
    
  }
  o <- merge(o, snp.dt, by="variant.id")
  return(o)  
                 }

  #### SAVE O
  output_file = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_omnibus_ALAN_MAY22023/"
  save(o.mods,
       file = paste(output_file,
                    "GLM_out.",
                    jobId,
                    ".",
                    "omnibus",
                    ".",
                    paste(wins.i$chr,wins.i$start,wins.i$end, sep = "_"),
                    ".Rdata", sep = ""))

  ####
  ####
  ####