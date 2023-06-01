# ijob -A berglandlab_standard -c20 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### get files
  fl <- list.files("/project/berglandlab/jcbnunez/Shared_w_Alan/GLMs_trialruns_Jun1_2023/", full.names=T)

  o <- foreach(fl.i=fl)%dopar%{
    # fl.i <- fl[1000]
    load(fl.i)
    message(paste(which(fl.i==fl), length(fl), sep=" / "))
    return(o.mods)

  }

  o <- rbindlist(o)

  table(o$chr)

### aggregate
  o.q <- o[,list(p_lrt=p_lrt, q=p.adjust(p_lrt, "fdr")), list(perm, model_features)]


  o.ag <- o.q[,list(nSig=mean(q<.05, na.rm=T)), list(perm, model_features)]

  o.ag.ag <- o.ag[,list(en=(nSig[perm==0])/median(nSig[perm!=0], na.rm=T)), list(model_features)]
  o.ag.ag
