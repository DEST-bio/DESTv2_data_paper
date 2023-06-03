# ijob -A berglandlab_standard -c4 -p dev --mem=10G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### get files
  fl <- list.files("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023", full.names=T)

### missing jobs
  fls <- gsub("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023/", "", fl)
  fls <- tstrsplit(fls, "\\.")[[2]]
  fls <- as.numeric(fls)
  c(1:9060)[!c(1:9060)%in%fls[order(fls)]]
  
    o <- foreach(fl.i=fl)%dopar%{
      # fl.i <- fl[100]
      load(fl.i)
      message(paste(which(fl.i==fl), length(fl), sep=" / "))
      return(o)

    }

  o <- rbindlist(o)

  table(o$chr)

### aggregate
  o.q <- o[,list(p_lrt=p_lrt, q=p.adjust(p_lrt, "fdr")), list(perm, model_features)]


  o.ag <- o[,list(nSig=mean(p_lrt<.005, na.rm=T)), list(perm, model_features)]

  o.ag.ag <- o.ag[,list(en=(nSig[perm==0])/median(nSig[perm!=0], na.rm=T)), list(model_features)]
  o.ag.ag
