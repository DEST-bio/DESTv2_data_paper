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

  #jobId=10

  permId <- jobId -1

### import all subpools
  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/", "contrast_perm_subpool_", full.name=T)

  c2p <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    load(fl.i)
    return(c2.dt.ag)
  }
  c2p <- rbindlist(c2p)
  setkey(c2p, perm)
  m.perm <- c2p[J(permId)]
  rm(c2p); gc()

### make window definitions
  win.bp <- 1e5
  step.bp <- 5e4

  setkey(m, "chr")

  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind",
                  .errorhandling="remove")%dopar%{

                    tmp <- m[J(chr.i)]

                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }

  wins[,i:=1:dim(wins)[1]]

  dim(wins)

### tack in observed data to get heterozygosity estimates
  load("~/dest2_glm_baypass_annotation.Rdata") ### made by `DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/collectBayPass.R`

  setkey(m, chr, pos)
  setkey(m.perm, chr, pos)

  m.perm <- merge(m.perm, m[,c("chr", "pos", "af"), with=F])
  m.perm[,cont.p:=10^(-neglog10P_median)]

### run windows
  m.perm[,C2.Z:=qnorm(cont.p, 0, 1)]
  m.perm[,het:=2*af*(1-af)]


  setkey(m.perm, chr, pos)
    win.out <- foreach(win.i=1:dim(wins)[1], .errorhandling = "remove", .combine = "rbind"  )%dopar%{
      # win.i <- 888
      message(paste(win.i, dim(wins)[1], sep=" / "))


      win.tmp <- m.perm[J(data.table(chr=wins[win.i]$chr,
                                      pos=wins[win.i]$start:wins[win.i]$end,
                                      key="chr,pos")), nomatch=0]


      win.out <- win.tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos),
                      win=win.i,
                      C2.wZa=sum(het*C2.Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                      C2.wZa.p=pnorm(sum(het*C2.Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                      nSNPs = .N,
                      perm=unique(perm)),]
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
  save(win.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/window_perms/C2_perm_", permId, ".Rdata", sep=""))
