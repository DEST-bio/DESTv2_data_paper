# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)

  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")

### get files
  fl <- list.files("/scratch/aob2x/dest_fstoutput/", "dest_fst_fet_window", full.name=T)

### collect completed jobs
  message("collect")
  o <- foreach(fl.i=fl, .errorhandling="pass")%dopar%{
    # fl.i <- fl[1]
    message(paste(which(fl.i==fl), length(fl), sep=" / "))

    load(fl.i)

    out[[1]]

  }
  o <- rbindlist(o)
  o_fst <- o

### outflank for each locyear
  of <- foreach(u=unique(o$locyear))%dopar%{
    #u <- unique(o$locyear)[1]
    tmp <- o[locyear==u]
    tmp[fet_p>1, fet_p:=.9999999]
    tmp[,fet_q:=qvalue(fet_p)$qvalue]
    tmp
    #OutFLANK(tmp, NumberOfSamples=2)$results
  }
  of <- rbindlist(of)
  of.ag <- of[,list(nSig_outFL=sum(pvaluesRightTail<1e-4)), list(variant.id)]
  of.ag[nSig_outFL>4]
  of[variant.id==2353046 & fet_q<.05]

  of[pvaluesRightTail<1e-5]

  of.ag <- of[,list(nSig_outFL=sum(fet_q<.1), st=mean(or[fet_q<.1]>1)), list(variant.id)]

  (table(of.ag[nSig_outFL==2]$st))
  binom.test(1527+1308, 5873, .5)


  load(file="/standard/vol186/bergland-lab/DEST_v2/GLMER_output_EuropeSeasonality.Rdata")
  o[,rnp:=rank(p_lrt)/(1+length(p_lrt))]

  o <- merge(o, o_fst, by="variant.id")


  o[rnp<.05, list(nSig=mean(fet_p<.05), or=mean(or[fet_p<.05]>1)), list(variant.id)]
