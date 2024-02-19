# ijob -A berglandlab_standard -c20 -p standard --mem=15G
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

### load data
  permId<- jobId - 1

  load(paste("/scratch/aob2x/DEST2_analysis/seasonality/perms/glmer_dest_perm", permId, ".Rdata", sep=""))
  m <- glmer.dest
  rm(glmer.dest)

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

### generate RNP values
  m[,rnp:=rank(p_lrt)/(1+length(p_lrt))]

### run windows
  setkey(m, chr, pos)

  win.out <- foreach(win.i=1:dim(wins)[1], .errorhandling = "remove", .combine = "rbind"  )%dopar%{
    # win.i <- 888
    message(paste(win.i, dim(wins)[1], sep=" / "))


    win.tmp <- m[J(data.table(chr=wins[win.i]$chr,
                                    pos=wins[win.i]$start:wins[win.i]$end,
                                    key="chr,pos")), nomatch=0]

    win.tmp <- win.tmp[!is.na(p_lrt)][p_lrt>0 & p_lrt<1]
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt, 0, 1)]

    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]

    win.tmp[,het:=2*af*(1-af)]

    pr.i =c(0.05)

    win.out <- win.tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos),
                  win=win.i, pr=pr.i,
                  glm.rnp.pr=c(mean(rnp<=pr.i)),
                  glm.rnp.binom.p=c(binom.test(sum(rnp<=pr.i),
                                           length(rnp), pr.i, alternative="greater")$p.value),
                  wZa=sum(het*Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                  wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                  min.p.lrt=min(p_lrt),
                  min.rnp=min(rnp),
                  nSNPs = .N,
                  perm=unique(win.tmp$perm)),]

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


  save(win.out, file=paste("/scratch/aob2x/DEST2_analysis/seasonality/window_perms/glm_perm_", permId, ".Rdata", sep=""))
