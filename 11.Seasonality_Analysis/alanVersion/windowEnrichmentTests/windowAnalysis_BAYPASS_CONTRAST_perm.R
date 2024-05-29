# ijob -A berglandlab_standard -c20 -p standard --mem=150G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(dplyr)

### get input variable
  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])

  #jobId=1

  permId <- jobId -1

### import all subpools
  # fl <- list.files("/scratch/aob2x/dest2_baypass/contrast_perms/", "contrast_perm_subpool_", full.name=T) ### this was totally random shuffling
  fl <- list.files("/scratch/aob2x/dest2_baypass/contrast_perms_v3", "contrast_perm_subpool_", full.name=T)
  fl <- fl[!grepl("RAW", fl)]
  c2p <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    load(fl.i)
    return(c2.dt.ag)
  }
  c2p <- rbindlist(c2p)
  setkey(c2p, perm)
  m.perm <- c2p
  #rm(c2p); gc()

### make window definitions
  win.bp <- 1e5
  step.bp <- 5e4

  setkey(m.perm, "chr")

  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind",
                  .errorhandling="remove")%dopar%{

                    tmp <- m.perm[J(chr.i)]

                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }

  wins[,i:=1:dim(wins)[1]]

  dim(wins)

### tack in observed data to get heterozygosity estimates
  load(file="~/dest2_glm_baypass_annotation_pod.podOutputToo_v3.Rdata")
  setnames(m, "C2_std_median", "C2_std_median_orig")
  setkey(m, chr, pos)
  setkey(m.perm, chr, pos)

  m.perm <- merge(m.perm, m[,c("chr", "pos", "af", "C2_std_median_orig"), with=F])
  m.perm[,cont.p:=10^(-neglog10P_median)]

  m.perm.ag <- m.perm[,list(pr=round(mean(round(C2_std_median, 4)<=round(C2_std_median_orig, 4)), 3)), list(perm)]
  m.perm.ag
  summary(m.perm.ag[perm!=0]$pr)

### load in ecdfs from PODs
  load(file="~/baypass_pod_ecdf.Rdata")
  m.perm[,C2.pod.p.orig:=1-cont.ecdf(m.perm$C2_std_median)]

  setkey(m.perm, perm)
  m.perm.ecdf <- foreach(i=0:100)%dopar%{
    #i <- 0
    tmp <- m.perm[J(i)]
    tmp.ecdf <- ecdf.list[[i+1]]
    tmp[, C2.pod.p:=1-tmp.ecdf(C2_std_median)]
  }
  m.perm.ecdf <- rbindlist(m.perm.ecdf)





### run windows
  m.perm[,C2.Z:=qnorm(cont.p, 0, 1)]
  m.perm[,C2.pod.Z:=qnorm(C2.pod.p, 0, 1)]

  m.perm[,het:=2*af*(1-af)]


  setkey(m.perm, chr, pos)

    win.out <- foreach(win.i=1:dim(wins)[1], .errorhandling = "remove", .combine = "rbind"  )%dopar%{
      # win.i <- 1900
      message(paste(win.i, dim(wins)[1], sep=" / "))


      win.tmp <- m.perm[J(data.table(chr=wins[win.i]$chr,
                                      pos=wins[win.i]$start:wins[win.i]$end,
                                      key="chr,pos")), nomatch=0]


      win.out <- win.tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos),
                               win=win.i,
                               C2.wZa=sum(het*C2.Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                               C2.wZa.p=pnorm(sum(het*C2.Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                               C2.wZa.pod=sum(het*C2.pod.Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                               C2.wZa.pod.p=pnorm(sum(het*C2.pod.Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                               nSNPs = .N),
                      list(perm=perm)]

        win.out[,invName:=case_when(
              chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
              chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
              chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
              chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
              chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
              chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
              TRUE ~ "noInv")]
        win.out


    }

### save
  win.out[,list(thr=quantile(C2.wZa.p, .01)), list(perm)][perm!=0][,list(thr=min(thr))]

  save(win.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/window_perms_podv3/C2_perm_", permId, ".Rdata", sep=""))




### check
  load("~/dest2_glm_baypass_annotation_pod.podOutputToo_v2.Rdata")
  setkey(m, chr, pos)
  setkey(m.perm, chr, pos)
  mm <- merge(m, m.perm)
  table(mm$C2_std_median.x==mm$C2_std_median.y)
