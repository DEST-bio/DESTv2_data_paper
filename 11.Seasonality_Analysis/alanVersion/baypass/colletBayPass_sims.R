# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)


### iterate
  fl <- list.files("/standard/vol186/bergland-lab/Gio/subsims", "xtx", full.name=T)
  xtx <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    file <- tstrsplit(fl.i, "/")%>%last()
    subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
    rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()
    simNum <- tstrsplit(file, "_")[[4]]%>%as.numeric()

    tmp <- fread(fl.i)
    setnames(tmp, "log10(1/pval)", "negLog10P")
    tmp[,subpool:=subPool]
    tmp[,rep:=rep]
    tmp[,simNum:=simNum]
    tmp
  }
  xtx <- rbindlist(xtx)
  xtx[,XtXst:=as.numeric(as.character(XtXst))]
  xtx.ag <- xtx[,list(XtXst_median=median(XtXst), M_P=median(M_P)),
                list(MRK, subpool, simNum)]
  xtx.ag[order(subpool)]
  xtx.ag[,list(q95=quantile(XtXst_median, .95, na.rm=T)), list(simNum)]

  save(xtx.ag, file="/scratch/aob2x/dest2_baypass_sims.Rdata")

### contrast
  fl <- list.files("/standard/vol186/bergland-lab/Gio/subsims", "contrast", full.name=T)
  contrast <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    file <- tstrsplit(fl.i, "/")%>%last()
    subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
    rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()
    simNum <- tstrsplit(file, "_")[[4]]%>%as.numeric()

    tmp <- fread(fl.i)
    setnames(tmp, "log10(1/pval)", "negLog10P")
    tmp[,subpool:=subPool]
    tmp[,rep:=rep]
    tmp[,simNum:=simNum]
    tmp
  }
  contrast <- rbindlist(contrast)

  contrast.ag <- contrast[,list(C2_std_median=median(C2_std)),
                list(MRK, subpool, simNum)]
  contrast.ag[,list(q95=quantile(C2_std_median, .999, na.rm=T)), list(simNum)]

  save(xtx.ag, file="/scratch/aob2x/dest2_baypass_sims.Rdata")
