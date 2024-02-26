# ijob -A berglandlab_standard -c20 -p standard --mem=15G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(dplyr)

### read data
  fl <- list.files("/scratch/aob2x/DEST2_analysis/seasonality/window_perms/", "glm_perm_",full.name=T)

  win.all <- foreach(fl.i=fl, .combine="rbind")%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    win.out
  }

  win.all.ag2 <- win.all[,list(wZa.p=wZa.p[perm==0],
                               perm.lci=quantile(wZa.p[perm!=0], .0),
                               perm.uci=quantile(wZa.p[perm!=0], 1),
                               perm.med=quantile(wZa.p[perm!=0], .5),
                               .N, noNA_n=length(!is.na(wZa.p)), chr=unique(chr), pos_min=median(pos_min), pos_max=median(pos_max), invName=invName[1]),
                          list(win)]
  win.all.ag2[,wZa.q:=p.adjust(exp(wZa.p), "bonferroni")]
  win.all.ag2[wZa.q<.05][wZa.p<perm.lci]

### save
    save(win.all, win.all.ag2, file="~/destv2_seasonality_perm.Rdata")

### local
  scp aob2x@rivanna.hpc.virginia.edu:~/destv2_seasonality_perm.Rdata ~/destv2_seasonality_perm.Rdata

  library(data.table)
  library(ggplot2)
  load("~/destv2_seasonality_perm.Rdata")
  win.all.ag2[,pos_mean:=pos_min/2 + pos_max/2]


    win.all.ag <- win.all[,list(min.wZa.p=min(wZa.p)), list(perm)]
    win.all.ag
    win.all.ag2[order(wZa.p)]



  ggplot(data=win.all.ag2, aes(x=pos_mean, y=-wZa.p)) +
  geom_line() + facet_grid(~chr) +
  geom_point(data=win.all.ag2[wZa.q<.05][wZa.p<perm.lci], aes(x=pos_mean, y=-wZa.p), color="red") +
  geom_hline(yintercept=-quantile(win.all.ag[perm!=0]$min.wZa.p, .05), color="red")

  quantile(win.all.ag[perm!=0]$min.wZa.p, .05)



  win.all.ag2 <- win.all[,list(wZa.p=wZa.p[perm==0], perm.thr=quantile(wZa.p[perm!=0], .05)), list(chr, pos_mean, invName)]
  table(win.all.ag2[wZa.p<perm.thr]$invName)
