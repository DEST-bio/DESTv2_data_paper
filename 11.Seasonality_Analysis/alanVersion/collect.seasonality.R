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

### jobs
  job.dt <- expand.grid(pops=c("NoCore20_seas", "NoCore20_NoProblems_seas"), mf=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))
  jobId=1:dim(job.dt)[1]

  # job.dt <- expand.grid(pops=c("Core20_bad"), mf=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))
  # jobId <- 1


### get files
  fns <- paste("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023", job.dt$pops[jobId], job.dt=job.dt$mf[jobId], sep="/")
  fl <- unlist(sapply(fns, list.files, full.names=T))
  length(fl)

### missing jobs
  # fls <- gsub("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023/", "", fl)
  # fls <- tstrsplit(fls, "\\.")[[2]]
  # fls <- as.numeric(fls)
  # # paste(c(1:9060)[!c(1:9060)%in%fls[order(fls)]], collapse=",")

### collect completed jobs
  o <- foreach(fl.i=fl, .errorhandling="stop")%dopar%{
    # fl.i <- fl[2]

    load(fl.i)
    message(paste(which(fl.i==fl), length(fl), sep=" / "))
    #return(o[pops==job.dt$pops[jobId]][model_features==job.dt$mf[jobId]])
    oo

  }
  o <- rbindlist(o, fill=T)

  table(o$chr, o$model_features, o$pops)
  table(is.na(o$p_lrt))
  o[,list(zero=sum( p_lrt==0, na.rm=T)), list(perm)]
  table(o$p_lrt==1)
  table(o[nObs>19]$p_lrt<.05, o[nObs>19]$perm, o[nObs>19]$model_features)
  table(o[nObs>90][model_features=="LocBinomial"]$p_lrt<.05, o[nObs>90][model_features=="LocBinomial"]$perm, o[nObs>90][model_features=="LocBinomial"]$chr)

  table(o$nFixed)

### save jobs based on model type
  setkey(o, model_features, pops)


  oo <- foreach(mf=unique(o$model_features), .combine="rbind")%do%{
    foreach(p=unique(o$pops), .combine="rbind")%do%{
      # mf="LocRan"; p="NoCore20_NoProblems_seas"
      message(paste("saving: ", mf, p, sep=" / "))
      mod.out <- o[J(data.table(model_features=mf, pops=p, key="model_features,pops"))]

      #save(mod.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/", p, "_", mf, ".Rdata", sep=""))
      setkey(mod.out, perm)

      #mod.out[,list(N=make_bins(p_lrt, size=.001, ret="my_sum"), thr=make_bins(p_lrt, size=.001, ret="my")), list(perm)]


      o.temp <- foreach(p.i=unique(mod.out$perm), .combine="rbind")%do%{
        foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
          # p.i <- 0
          message(p.i)

          #o.sig <- make_bins(x=mod.out[J(p.i)]$p_lrt, size=.001)
          #data.table(nSig=o.sig[[1]], thr=o.sig[[2]][-1], perm=p.i, model_features=mf, pops=p)
          grid <- 0.001
          my_seq = data.table(min_p=seq(from=0, to=1-grid, by=grid), max_p=seq(from=grid, to=1, by=grid))

          tmp <- mod.out[nObs>80][J(p.i)][!is.na(p_lrt)][][chr==chr.i][nFixed==0][af>.1 & af<.9][my_seq, .(N = .N), on = .(p_lrt > min_p, p_lrt <= max_p), by = .EACHI]
          setnames(tmp, c("min_p", "max_p", "N"))
          tmp[,perm:=p.i]
          tmp[,pops:=p]
          tmp[,model_features:=mf]
          tmp[,chr:=chr.i]
          tmp
        }
      }
      return(o.temp)

    }
  }
  #oo[max_p==.001][perm<=2][order(model_features)][chr=="2R"]

  save(oo, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.", job.dt$pops[2], ".Rdata", sep=""))



### sliding window aggregation
  o.ag <- o[,list(nSig=sum(p_lrt<1e-4)), list(model_features, chr, pos, perm=perm!=0)]
  o.ag[perm==T & nSig>8]
  o.ag[perm==F & nSig>0]
