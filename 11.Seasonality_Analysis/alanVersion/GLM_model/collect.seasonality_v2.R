# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
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
  o.perm.quan <- o[,list(q999=quantile(-log10(p_lrt), .999)), list(perm)]


  o <- o[perm==0][nObs>100][af>.05 & af<.95][neff>20][,list(chr,pos,invName,variant.id,p_lrt)]
  save(o, file="/standard/vol186/bergland-lab/DEST_v2/GLMER_output_EuropeSeasonality.Rdata")

### some summaries
  o.q <- o[,list(q=p.adjust(p_lrt, "fdr"), rnp=rank(p_lrt)/(1+length(p_lrt)), chr, pos), list(perm)]
  setkey(o, chr, pos, perm)
  setkey(o.q, chr, pos, perm)
  o <- merge(o, o.q)
  #o.ag <- o[,list(nSig=sum( p_lrt<.0005, na.rm=T), nSingular=sum(singular=="TRUE_TRUE"), p1=sum(p_lrt>.99, na.rm=T)), list(perm, model_features,chr,inv=invName!="noInv")]
  o.ag <- o[,list(nSig_p=sum( p_lrt<.0005, na.rm=T),
                  nRNP=sum(rnp<.05),
                  nSig_q=sum(q<.05, na.rm=T)), list(perm, model_features,chr,inv=invName!="noInv")]

  o.ag[,list(pr_p=mean(nSig_p[perm==0]>nSig_p[perm!=0]),
            pr_rnp=mean(nRNP[perm==0]>nRNP[perm!=0]),
            pr_q=mean(nSig_q[perm==0]>nSig_q[perm!=0])), list(chr,inv)]

### saving the
  o.small <- o[perm<=3]
  save(o.small, file="/scratch/aob2x/o_small.Rdata")

#### enrichment
### save jobs based on model type
  setkey(o, model_features, pops)

  oo <- foreach(mf=unique(o$model_features), .combine="rbind")%do%{
    foreach(p=unique(o$pops), .combine="rbind")%do%{
      # mf=unique(o$model_features); p=unique(o$pops); p.i=0; chr.i="2L"; inv.i=T
      message(paste("saving: ", mf, p, sep=" / "))
      mod.out <- o[J(data.table(model_features=mf, pops=p, key="model_features,pops"))]

    #  message(head(mod.out))

      # folder_out = paste(folder, "/compiled/glm_output/", sep = "")
      # save(mod.out, file=paste(folder_out, p, "_", mf, ".Rdata", sep=""))
       setkey(mod.out, perm)

      #mod.out[,list(N=make_bins(p_lrt, size=.001, ret="my_sum"), thr=make_bins(p_lrt, size=.001, ret="my")), list(perm)]

      nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T)
      o.temp <- foreach(p.i=unique(mod.out$perm), .combine="rbind")%do%{
        foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
          foreach(inv.i=c(TRUE, FALSE), .combine="rbind")%do%{
            # p.i <- 0; chr.i<-"2L"; inv.i<-T
            message(paste(p.i, chr.i, inv.i, sep=" / "))

            ### tabluate
              grid <- 0.001
              my_seq = data.table(min_p=seq(from=0, to=1-grid, by=grid), max_p=seq(from=grid, to=1, by=grid))

              tmp <- mod.out[(invName!="noInv")==inv.i][nObs>=nObs.thr][J(p.i)][!is.na(p_lrt)][chr==chr.i][nFixed==0][af>.05 & af<.95][my_seq, .(N = .N), on = .(p_lrt > min_p, p_lrt <= max_p), by = .EACHI]

            ### format and return
              setnames(tmp, c("min_p", "max_p", "N"))
              tmp[,perm:=p.i]
              tmp[,pops:=p]
              tmp[,model_features:=mf]
              tmp[,chr:=chr.i]
              tmp[,inv:=inv.i]
              tmp
            }
          }
      }
      # o.temp[max_p==.005][chr=="2L"]
      return(o.temp)

    }
  }
  #oo[max_p==.001][perm<=2][order(model_features)][chr=="2L"][pops=="NoCore20_NoProblems_Steep_Neg_seas"][chr=="2L"][inv==T]
  head(oo)


  save(oo, file=paste(paste(folder, "/compiled/", sep = ""),"/enrichment.", p, "_", mf, ".Rdata", sep=""))
