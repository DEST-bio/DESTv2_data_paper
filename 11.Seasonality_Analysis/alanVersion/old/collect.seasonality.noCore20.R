# ijob -A berglandlab_standard -c5 -p dev --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### get jobs
  args = commandArgs(trailingOnly=TRUE)

  jobId=as.numeric(args[1])
  #jobId=1

### jobs
  job.dt <- expand.grid(pops="Core20_seas", mf=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))


### get files
#  fl <- list.files(paste("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023", job.dt$pops[jobId], job.dt=job.dt$mf[jobId], sep="/"), full.names=T)
  fl <- list.files("/project/berglandlab/jcbnunez/Shared_w_Alan/GLMs_trialruns_Jun1_2023", full.names=F)

### missing jobs
  # fls <- gsub("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023/", "", fl)
  # fls <- tstrsplit(fls, "\\.")[[2]]
  # fls <- as.numeric(fls)
  # # paste(c(1:9060)[!c(1:9060)%in%fls[order(fls)]], collapse=",")

### collect completed jobs
  o <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[100]
    load(fl.i)
    message(paste(which(fl.i==fl), length(fl), sep=" / "))
    return(o.mods[perm<=10][!is.na(p_lrt) & p_lrt!=0])

  }
  o <- rbindlist(o, fill=T)

  table(o$chr)

### save jobs based on model type
  setnames(o, "model", "pops")
  setkey(o, model_features, pops)


  oo <- foreach(mf=unique(o$model_features), .combine="rbind")%do%{
    foreach(p=unique(o$pops), .combine="rbind")%do%{
      # mf="No_Phylo"; p="NoCore20_seas"
      message(paste("saving: ", mf, p, sep=" / "))
      mod.out <- o[J(data.table(model_features=mf, pops=p, key="model_features,pops"))]

      #save(mod.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled", p, "_", mf, ".Rdata", sep=""))
      setkey(mod.out, perm)

      #mod.out[,list(N=make_bins(p_lrt, size=.001, ret="my_sum"), thr=make_bins(p_lrt, size=.001, ret="my")), list(perm)]


      o.temp <- foreach(p.i=unique(mod.out$perm), .combine="rbind")%do%{
        # p.i <- 0
        message(p.i)

        #o.sig <- make_bins(x=mod.out[J(p.i)]$p_lrt, size=.001)
        #data.table(nSig=o.sig[[1]], thr=o.sig[[2]][-1], perm=p.i, model_features=mf, pops=p)
        my_seq = data.table(min_p=seq(from=0, to=1-.001, by=.001), max_p=seq(from=0.001, to=1, by=.001))
        tmp <- mod.out[J(p.i)][my_seq, .(N = .N), on = .(p_lrt > min_p, p_lrt < max_p), by = .EACHI]
        setnames(tmp, c("min_p", "max_p", "N"))
        tmp[,perm:=p.i]
        tmp[,pops:=p]
        tmp[,model_features:=mf]

        tmp
      }
      return(o.temp)
    }
  }
  oo[thr==.001][model_features=="LocBinomial"]

  save(oo, file="/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.NoCore20_seas.Rdata")
