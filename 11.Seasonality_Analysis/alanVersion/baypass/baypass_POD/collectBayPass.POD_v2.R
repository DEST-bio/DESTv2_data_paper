
# ijob -A berglandlab -c5 -p largemem --mem=250G
# ijob -A berglandlab_standard -c20 -p standard --mem=250G

### module load gcc/11.4; module load openmpi/4.1.4; module load R/4.3.1; R
### libraries
  .libPaths(c("/project/berglandlab/Rlibs_4.3.1/")); .libPaths()

  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(dplyr)
  library(SeqArray)

#########
### XtX
########


  ### get files
    fl <- list.files("/scratch/aob2x/dest2_baypass/pods_v2/", "xtx", full.name=T)

  ### iterate
    xtx <- foreach(fl.i=fl)%dopar%{
      message(fl.i)
      # fl.i <- fl[1]
      file <- tstrsplit(fl.i, "/")%>%last()
      subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
      rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()

      tmp <- fread(fl.i)
      setnames(tmp, "log10(1/pval)", "negLog10P")
      tmp[,subpool:=subPool]
      tmp[,rep:=rep]
      tmp
    }
    xtx <- rbindlist(xtx)
    xtx[,XtXst:=as.numeric(as.character(XtXst))]

  ## merge
    #setkey(snp.dt, MRK, subpool)
    #setkey(xtx, MRK, subpool)
    #xtx <- merge(xtx, snp.dt)
    #save(xtx, file="/scratch/aob2x/xtx_dest2.pod.raw.Rdata")

  ### average
    xtx[negLog10P==Inf, negLog10P:=max(xtx[negLog10P!=Inf]$negLog10P)]

    xtx.ag <- xtx[,list(XtXst_median=median(XtXst), XtXst_mean=mean(XtXst),
                        xtx_neglogp_median=median(negLog10P), xtx_neglogp_mean=mean(negLog10P), af=mean(M_P)), list(MRK, subpool)]

##########
### contrast
##########
  ### get files
    fl <- list.files("/scratch/aob2x/dest2_baypass/pods_v2/", "contrast", full.name=T)

    cont <- foreach(fl.i=fl)%dopar%{
      message(fl.i)
      # fl.i <- fl[1]
      file <- tstrsplit(fl.i, "/")%>%last()
      subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
      rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()

      tmp <- fread(fl.i)
      setnames(tmp, "log10(1/pval)", "negLog10P")
      tmp[,subpool:=subPool]
      tmp[,rep:=rep]
      tmp
    }
    cont <- rbindlist(cont)

  ### average
    cont.ag <- cont[,list(C2_std_median=median(C2_std), C2_std_mean=mean(C2_std),
                          C2_neglogp_median=median(negLog10P), C2_neglogp_mean=mean(negLog10P)), list(MRK, subpool)]

### load in empirical data
  load(file="~/dest2_glm_baypass_annotation.Rdata")


### quick check: calculate empirical p-values based on permutations
  summary(xtx.ag$XtXst_median)
  summary(m$XtXst_median)

  summary(cont.ag$C2_std_median)
  summary(m$C2_std_median)

### calcualte POD based thresholds
  xtx.ecdf <- ecdf(xtx.ag$XtXst_median)
  m[,XtXst.pod.p :=  1 - xtx.ecdf(XtXst_median)]
  m[XtXst.pod.p==0, XtXst.pod.p:=1/dim(m)[1]]


  cont.ecdf <- ecdf(cont.ag$C2_std_median)
  m[,cont.pod.p :=  1 - cont.ecdf(C2_std_median)]
  m[cont.pod.p==0, cont.pod.p:=1/dim(m)[1]]

  m[which.min(xtx.p)]
  m[which.min(XtXst.pod.p)]

  summary(m[min(XtXst.pod.p)==XtXst.pod.p])
  summary(m[min(cont.pod.p)==cont.pod.p])

### thresholds
  thrs <- c(.95, .99, .995, .999, .9995, .9999)
  thrs.ag <- cbind(
  xtx.ag[,list(XtXst_thr=quantile(XtXst_median, thrs, na.rm=T), thr=thrs)],
  cont.ag[,list(C2_thr=quantile(C2_std_median,  thrs, na.rm=T))])

### save
  save(m, thrs.ag, xtx.ag, cont.ag, file="~/dest2_glm_baypass_annotation_pod.podOutpuToo.Rdata")

summary(lm(-log10(cont.p)~I(-log10(cont.pod.p)), data=m)
