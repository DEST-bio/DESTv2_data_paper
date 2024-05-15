# ijob -A berglandlab_standard -c10 -p largemem --mem=150G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(dplyr)

### import all subpools
  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/", "contrast_perm_subpool_", full.name=T)

  c2p <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    load(fl.i)
    return(c2.dt.ag)
  }
  c2p <- rbindlist(c2p)

### outlier test
  c2p.outlier <- c2p[,list(maxc=max(C2_std_median), nsig=mean(C2_std_median>3)), list(perm)]
  c2p.outlier[order(nsig)]






### li
  library(data.table)
  library(foreach)
  library(ggplot2)

### data
  load("~/dest2_glm_baypass_annotation.Rdata")
  prop.table(table(m$C2_std_median>3))

  hist(m$C2_neglogp_median)
  m[,cont.p2:=1-pchisq(C2_std_mean, 1)]
  hist(m$cont.p2)



  dat <- fread("/Users/alanbergland/standard/destsubpool_21_2_summary_contrast.out")
  dat[,cont.p2:=1-pchisq(C2_std, 1)]
  hist(dat$cont.p2)
