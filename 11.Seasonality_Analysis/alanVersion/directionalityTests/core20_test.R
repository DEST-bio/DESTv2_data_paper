# ijob -A berglandlab_standard -c5 -p dev --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)

### wd
  setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_13_2023/compiled/glm_output")

### load in datasets
  ### focal
    #load("NoCore20_seas_Phylo_yearPop_Ran.Rdata")
    load("Core20_seas_Phylo_yearPop_Ran.Rdata")

    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    focal <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

  ### original core20
    core20.orig <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")
    core20.swap <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

    core20.orig[,set:="orig"]
    core20.swap[,set:="swap"]
    core20 <- rbind(core20.orig, core20.swap)

    setnames(core20, "chrom", "chr")

  ### R5 -> R6 DGRP conversion table
    liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
    liftover <- fread(liftover.fn)
    liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

  ### do liftover
    setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
    setkey(core20, dm3_chr, dm3_pos)
    setkey(liftover, dm3_chr, dm3_pos)

    core20 <- merge(core20, liftover)

    setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

  ### merge with VA_glm
    setkey(core20, chr, pos)
    setkey(focal, chr, pos)
    m <- merge(focal[perm==0], core20)

### basic test
  m$rnp.new <- rank(m$p_lrt)/(1+length(m$p_lrt))
  m$rnp.old <- rank(m$seas.p)/(1+length(m$seas.p))

  fisher.test(table(m$rnp.new<.05, m$rnp.old<.05))

  prop.test(table(sign(m[rnp.new<.01 & rnp.old<.01]$b_seas) == sign(m[rnp.new<.01 & rnp.old<.01]$seas.coef) ))







  ### clean up
    rm(mod.out)

### merge
  setkey(focal, chr, pos, perm)
  setkey(tester, chr, pos, perm)

  m <- merge(focal, tester[,c("chr", "pos", "perm", "p_lrt", "b_seas"), with=F])

### perms
  setkey(m, perm)
  foreach(i=0:10, .combine="rbind")%dopar%{
    tmp <- m[J(i)]
    tmp[,rnp.x:=rank(p_lrt.x)/(length(p_lrt.x)+1)]
    tmp[,rnp.y:=rank(p_lrt.y)/(length(p_lrt.y)+1)]

    #isher.test(tmp[chr=="2L"]$rnp.x<.01, tmp[chr=="2L"]$rnp.y<.01)

    tmp2 <- tmp[chr=="2L"][rnp.x<.05 & rnp.y<0.05]
    data.table(perm=i, prop=prop.test(rev(table(sign(tmp2$b_seas.y)==sign(tmp2$b_seas.x))))$estimate)
  }
