# ijob -A berglandlab_standard -c5 -p largemem --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)

### wd

### load in datasets
  ### new
    setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_13_2023/compiled/glm_output")

    load("NoCore20_NoProblems_Steep_Neg_seas_yearPop_Ran.Rdata")
    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    newr <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

  ### old
    setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/glm_output")

    load("NoCore20_NoProblems_NoFlip_seas_LocRan.Rdata")
    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    oldr <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

### merge
  setkey(newr, chr, pos, perm)
  setkey(oldr, chr, pos, perm)

  m <- merge(newr, oldr, all=T)

  table(m$chr, is.na(m$p_lrt.x))

  table(m[chr=="2L"]$p_lrt.x < 1e-3,  m[chr=="2L"]$p_lrt.y <1e-3, m[chr=="2L"]$perm)


###
  table(newr$p_lrt<1e-3, newr$perm, newr$chr)
  table(oldr$p_lrt<1e-4, oldr$perm, oldr$chr)







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
