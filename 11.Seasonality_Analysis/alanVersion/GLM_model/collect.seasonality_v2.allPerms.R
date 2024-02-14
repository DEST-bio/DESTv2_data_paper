# ijob -A berglandlab_standard -c5 -p largemem --mem=250G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)


### get jobs
  args = commandArgs(trailingOnly=TRUE)

  jobId=as.numeric(args[1])
  message(jobId)

### paths
  folder <- "/scratch/aob2x/DEST2_analysis/seasonality/glm_test_SEPT_29_2023/"
  pops <- "NoCore20_seas_europe"
  method <- "yearPop_Ran"

### get files
  fns <- paste(folder, pops, method, sep="/")

  fl <- as.vector(unlist(sapply(fns, list.files, full.names=T)))
  length(fl)

### missing jobs
  # fls <- gsub("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_1_2023/", "", fl)
  # fls <- tstrsplit(fls, "\\.")[[2]]
  # fls <- as.numeric(fls)
  # # paste(c(1:9060)[!c(1:9060)%in%fls[order(fls)]], collapse=",")

### collect completed jobs
  message("collect")
  o <- foreach(fl.i=fl, .errorhandling="pass")%dopar%{
    # fl.i <- fl[1]
    message(paste(which(fl.i==fl), length(fl), sep=" / "))

    load(fl.i)

    if(grepl("Core20_bad", fl.i)){
      oo <- o[model_features==tstrsplit(fl.i, "/")[[8]]]
    }
    oo[,invName:=case_when(
          chr=="2L" & pos >	2225744	 & pos < 13154180	 ~ "2Lt",
          chr=="2R" & pos >	15391154 & pos < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos >	11750567 & pos < 	26140370 ~ 	"3RK",
          chr=="3R" & pos >	21406917 & pos < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos >	16432209 & pos < 	24744010 ~ 	"3RP",
          chr=="3L" & pos >	3173046	 & pos < 16308841	 ~ "3LP",
          TRUE ~ "noInv")]
    return(oo)
    #oo[perm==0]

  }
  o <- rbindlist(o, fill=T)
  o <- o[nObs>100][af>.05 & af<.95][neff>20]
  #o.perm.quan <- o[,list(q999=quantile(-log10(p_lrt), .999)), list(perm)]


  #o <- o[nObs>100][af>.05 & af<.95][neff>20][,list(chr,pos,invName,variant.id,p_lrt, perm)]
  setkey(o, perm)

  foreach(i=0:100)%do%{
    # i<-0
    glmer.dest <- o[J(i)]
    save(glmer.dest, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/perms/glmer_dest_perm", i, ".Rdata", sep=""))
  }
