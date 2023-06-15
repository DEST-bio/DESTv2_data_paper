# ijob -A berglandlab_standard -c5 -p dev --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load data

### wd
  setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/glm_output")

### load in datasets
  ### focal
    load("NoCore20_NoProblems_NoFlip_seas_LocRan.Rdata")
    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    mod.out <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

### make window definitions
  win.bp <- 1e5
  step.bp <- 5e4

  setkey(mod.out, "chr")

  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind",
                  .errorhandling="remove")%dopar%{

                    tmp <- mod.out[J(chr.i)]

                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }

  wins[,i:=1:dim(wins)[1]]

  dim(wins)

### generate RNP values
  mod.out.rnp <- mod.out[,list(rnp=rank(p_lrt)/(length(p_lrt)+1), variant.id), list(model_features, perm)]

  setkey(mod.out.rnp, variant.id, perm, model_features)
  setkey(mod.out, variant.id, perm, model_features)

  mod.out <- merge(mod.out, mod.out.rnp)

### run windows

  setkey(mod.out, chr, pos)

  win.out <- foreach(win.i=1:dim(wins)[1], .errorhandling = "remove", .combine = "rbind"  )%dopar%{
    # win.i <- 736
    message(paste(win.i, dim(wins)[1], sep=" / "))


    win.tmp <- mod.out[J(data.table(chr=wins[win.i]$chr,
                                    pos=wins[win.i]$start:wins[win.i]$end,
                                    key="chr,pos")), nomatch=0]

    win.tmp <- win.tmp[!is.na(p_lrt)][p_lrt>0 & p_lrt<1]
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]

    win.tmp[,het:=2*af*(1-af)]

    pr.i <- c( 0.05)

    win.out <- win.tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos),
                  win=win.i, pr=pr.i,
                  rnp.pr=c(mean(rnp<=pr.i)),
                  rnp.binom.p=c(binom.test(sum(rnp<=pr.i),
                                           length(rnp), pr.i)$p.value),
                  wZa=sum(het*Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                  wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                  rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                  rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                  min.p.lrt=min(p_lrt),
                  min.rnp=min(rnp),
                  nSNPs = .N,
                  sum.rnp=sum(rnp<=pr.i)), list(model_features, perm)]

    win.out[, perm_type:=ifelse(perm==0, "real","permuted")]
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

  save(win.out, file="~/NoCore20_NoProblems_NoFlip_seas_LocRan.windows.Rdata")


  table(win.out$wZa.p<1e-10, win.out$perm)




### plot
  system("scp aob2x@rivanna.hpc.virginia.edu:~/NoCore20_NoProblems_NoFlip_seas_LocRan.windows.Rdata ~/.")

  load("~/NoCore20_NoProblems_NoFlip_seas_LocRan.windows.Rdata")

  ggplot(data=win.out[perm!=0]) +
    geom_line(aes(x=pos_mean, y=-log10(wZa.p), group=perm), alpha=.5) +
    geom_line(data=win.out[perm==0],aes(x=pos_mean, y=-log10(wZa.p)), color="red") +
    facet_grid(~chr)

  win.out.pa <- win.out[,list(wZa.pa=p.adjust(wZa.p, "bonferroni"), win), list(perm)]

  setkey(win.out, perm, win)
  setkey(win.out.pa, perm, win)
  win.out <- merge(win.out, win.out.pa)


  win.out.ag <- win.out[,list()]
