# ijob -A berglandlab_standard -c4 -p dev --mem=10G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(2)
  
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

  make_bins = function(x, size){
      x = na.omit(x)
      my_seq = seq(from=0, to=1, by=size)
      my_sum=vector()
      for (i in 1:(length(my_seq)-1) ){
          my_sum[i] = sum(x>=my_seq[i] & x<my_seq[i+1])
      }
      list(my_sum, my_seq)
  }


  oo <- foreach(mf=unique(o$model_features), .combine="rbind")%dopar%{
    foreach(p=unique(o$pops), .combine="rbind")%do%{
      # mf="Loc"; p="Core20_seas"
      message(paste("saving: ", mf, p, sep=" / "))
      mod.out <- o[J(data.table(model_features=mf, pops=p, key="model_features,pops"))]

      save(mod.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/", p, "_", mf, ".Rdata", sep=""))

      foreach(p.i=unique(mod.out$perm), .combine="rbind")%do%{
        # p.i <- 0
        o.sig <- make_bins(x=mod.out[perm==p.i]$p_lrt, size=.001)
        data.table(nSig=o.sig[[1]], thr=o.sig[[2]][-1], perm=p.i, model_features=mf, pops=p)
      }

    }
  }

  save(oo, file="/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata")
