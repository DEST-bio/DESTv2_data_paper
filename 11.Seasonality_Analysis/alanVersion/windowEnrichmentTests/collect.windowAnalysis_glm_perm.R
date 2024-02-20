# ijob -A berglandlab_standard -c20 -p standard --mem=15G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(dplyr)

### read data
  fl <- list.files("/scratch/aob2x/DEST2_analysis/seasonality/window_perms/", "glm_perm_",full.name=T)

  win.all <- foreach(fl.i=fl, .combine="rbind")%dopar%{
    load(fl.i)
    win.out
  }


### save
  save(win.all, file="~/destv2_seasonality_perm.Rdata")


scp aob2x@rivanna.hpc.virginia.edu:~/destv2_seasonality_perm.Rdata ~/destv2_seasonality_perm.Rdata

  library(data.table)
  library(ggplot2)
  load("~/destv2_seasonality_perm.Rdata")

  ggplot(data=win.all, aes(x=pos_mean, y=wZa.p, group=perm, color=as.factor(perm==0))) + geom_line() + facet_grid(~chr)


  win.all.ag <- win.all[,list(min.wZa.p=min(wZa.p)), list(perm)]
  win.all.ag
  win.all.ag[order(min.wZa.p)]

  quantile(win.all.ag[perm!=0]$min.wZa.p, .05)



  win.all.ag2 <- win.all[,list(wZa.p=wZa.p[perm==0], perm.thr=quantile(wZa.p[perm!=0], .05)), list(chr, pos_mean, invName)]
  table(win.all.ag2[wZa.p<perm.thr]$invName)
