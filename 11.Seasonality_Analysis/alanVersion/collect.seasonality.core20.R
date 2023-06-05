# ijob -A berglandlab_standard -c4 -p dev --mem=10G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### get files
  fl <- list.files("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023", full.names=T)

### missing jobs
  fls <- gsub("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023/", "", fl)
  fls <- tstrsplit(fls, "\\.")[[2]]
  fls <- as.numeric(fls)
  # paste(c(1:9060)[!c(1:9060)%in%fls[order(fls)]], collapse=",")

### collect completed jobs
  o <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[100]
    load(fl.i)
    message(paste(which(fl.i==fl), length(fl), sep=" / "))
    return(o[!is.na(p_lrt) & p_lrt!=0])

  }
  o <- rbindlist(o, fill=T)

  table(o$chr)

### save jobs based on model type
  setkey(o, model_features, pops)
  foreach(mf=unique(o$model_features))%dopar%{
    foreach(p=unique(o$pops))%do%{
      # mf="Loc"; p="Core20_seas"
      mod.out <- o[J(data.table(model_features=mf, pops=p, key="model_features,pops"))]
      save(mod.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/", p, "_", mf, ".Rdata", sep=""))
    }
  }

### aggregate
  pthr <- c(0, expand.grid(sapply(c(-6:-1), function(x) c(1:9)*10^x))[,1])

  o.sig <- foreach(x=1:(length(pthr)-1))%do%{
    # x <- pthr[10,1]
    message(x)
    o[,list(nSig=sum(p_lrt>pthr[x] & p_lrt<=pthr[x+1]),
            thr_min=pthr[x], thr_max=pthr[x+1],
            N=.N),
       list(perm, model_features, chr)]

  }
  o.sig <- rbindlist(o.sig)
  save(o.sig, file="/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas")
