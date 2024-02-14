# ijob -A berglandlab_standard -c5 -p dev --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### wd
  setwd("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_13_2023/compiled/glm_output")

### load in datasets
    load("Core20_seas_yearPop_Ran.Rdata")
    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    core20 <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

  ### old
    load("NoCore20_NoProblems_Steep_Neg_seas_yearPop_Ran.Rdata")
    nObs.thr <- quantile(mod.out$nObs, .25, na.rm=T); nObs.thr
    newSeas <- mod.out[nObs>=nObs.thr][!is.na(p_lrt)][nFixed==0][af>.05 & af<.95]

### merge
  setkey(core20, chr, pos, perm)
  setkey(newSeas, chr, pos, perm)
  m <- merge(core20, newSeas)

### load in original core20
  core20.orig <- fread("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")

  core20.orig[,set:="orig"]
  core20 <- rbind(core20.orig)

  setnames(core20, "chrom", "chr")

  ### R5 -> R6 DGRP conversion table
  liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
  liftover <- fread(liftover.fn)
  liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

  ### do liftover
  setnames(core20, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
  setkey(core20, dm3_chr, dm3_pos)
  setkey(liftover, dm3_chr, dm3_pos)

  core20 <- merge(core20, liftover)

  setnames(core20, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

  setkey(core20, chr, pos)
  setkey(m, chr, pos)
  mm <- merge(core20, m, all=T)

### save
  save(mm, file="~/NoCore20_NewCore20_OrigCore20.Rdata")

### load
# ijob -A berglandlab_standard -c5 -p dev --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  library(dplyr)
  registerDoMC(5)

  load(file="~/NoCore20_NewCore20_OrigCore20.Rdata")


### simple tests
  table(mm$seas.p<1e-3, mm$p_lrt.x<1e-3, mm$perm)
  table(mm$seas.p<1e-3, mm$p_lrt.y<1e-3, mm$perm)
  table(mm[!is.na(seas.p)]$p_lrt.x<.05, mm[!is.na(seas.p)]$p_lrt.y<.05, mm[!is.na(seas.p)]$perm)
  tab <- table(mm$p_lrt.x<.05, mm$p_lrt.y<.05, mm$perm)

### make window definitions
  win.bp <- 100000
  step.bp <- 50000

  setkey(mm, "chr")

  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind",
                  .errorhandling="remove")%dopar%{

                    tmp <- mm[J(chr.i)]

                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }

  wins[,i:=1:dim(wins)[1]]
  summary(wins$start - wins$end, na.rm=T)
  dim(wins)


### generate RNP values
  seas.rnp <- mm[!is.na(seas.p) & seas.p>0 & seas.p<1,    list(seas.p.rnp=rank(seas.p)/(length(seas.p)+1), variant.id.x), list(perm)]
  p_lrt.x.rnp <- mm[!is.na(p_lrt.x) & p_lrt.x>0 & p_lrt.x<1,list(p_lrt.x.rnp=rank(p_lrt.x)/(length(p_lrt.x)+1), variant.id.x), list(perm)]
  p_lrt.y.rnp <- mm[!is.na(p_lrt.y) & p_lrt.y>0 & p_lrt.y<1,list(p_lrt.y.rnp=rank(p_lrt.y)/(length(p_lrt.y)+1), variant.id.x), list(perm)]


  setkey(mm, variant.id.x, perm)
  setkey(seas.rnp, variant.id.x, perm)
  setkey(p_lrt.x.rnp, variant.id.x, perm)
  setkey(p_lrt.y.rnp, variant.id.x, perm)

  mm <- merge(mm, p_lrt.x.rnp)
  mm <- merge(mm, p_lrt.y.rnp)
  mm2 <- merge(mm, seas.rnp[!is.na(perm)])

  setkey(mm2, perm)

  tab <- table(mm2[J(0)]$seas < .05, mm2[J(0)]$p_lrt.x< .05); fisher.test(tab)
  tab <- table(mm2[J(0)]$seas < .05, mm2[J(0)]$p_lrt.y< .05); fisher.test(tab)
  tab <- table(mm2[J(0)]$seas.rnp < .05, mm2[J(0)]$p_lrt.x.rnp< .05); fisher.test(tab)
  tab <- table(mm2[J(0)]$seas.rnp < .05, mm2[J(0)]$p_lrt.y.rnp< .05); fisher.test(tab)
  tab <- table(mm2[J(0)]$p_lrt.x.rnp < .05, mm2[J(0)]$p_lrt.y.rnp< .05); fisher.test(tab)

### generate Z values
  mm2[,seas.p.Z:=qnorm(seas.p, 0, 1)]
  mm2[,p_lrt.x.Z:=qnorm(p_lrt.x, 0, 1)]
  mm2[,p_lrt.y.Z:=qnorm(p_lrt.y, 0, 1)]

  mm2[,het.core20:=2*af.x*(1-af.x)]
  mm2[,het.noCore20:=2*af.y*(1-af.y)]

### combo p-value
  mm2[,combo_p:=1-pchisq(-2*(log(p_lrt.x) + log(p_lrt.y)), 4)]
  mm2[,combo_p_lrt:=1-pchisq(-2*(log(p_lrt.x.rnp) + log(p_lrt.y.rnp)), 4)]
  table(mm2$perm, mm2$combo_p_lrt<1e-3)

  thr <- 1e-2
  table(sign(mm2[combo_p_lrt<thr]$b_seas.x) == sign(mm2[combo_p_lrt<thr]$b_seas.y), mm2[combo_p_lrt<thr]$perm)


  thr <- .01
  mm2[combo_p<1e-5, list(dir=mean(sign(b_seas.y)==sign(b_seas.x)), .N), list(perm)][order(dir)]
  d


  ggplot(data=mm2, aes(x=combo_p)) + geom_histogram(bins=100) + facet_grid(~perm)

### run windows

  setkey(mm2, chr, pos)

  win.out <- foreach(win.i=1:dim(wins)[1], .errorhandling = "remove", .combine = "rbind"  )%dopar%{
    # win.i <- 736
    message(paste(win.i, dim(wins)[1], sep=" / "))

    ### do calculations
      pr.i <- c( 0.05)
      foreach(term.i=c("seas.p", "p_lrt.x", "p_lrt.y"), .combine="rbind")%do%{
        # term.i <- "p_lrt.y"
        tmp <- mm2[J(data.table(chr=wins[win.i]$chr,
                                pos=wins[win.i]$start:wins[win.i]$end,
                                key="chr,pos")), nomatch=0]
        setnames(tmp, paste(term.i, "Z", sep="."), "Z")
        setnames(tmp, paste(term.i, "rnp", sep="."), "rnp")
        if(term.i!="p_lrt.y") {
          setnames(tmp, "het.core20", "het")
        } else if(term.i=="p_lrt.y") {
          setnames(tmp, "het.noCore20", "het")
        }

        win.out <- tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos), win=win.i, pr=pr.i,
                      rnp.pr=mean(rnp<=pr.i),
                      wZa=sum(het*Z, na.rm=T)/sqrt(sum(het^2, na.rm=T)),
                      nSNPs = .N),list(perm)]
        win.out[,term:=term.i]
      }
  }

  win.out[,wZa.p:=pnorm(wZa, lower.tail=T)]
  win.out[,rnp.p:=pbinom(rnp.pr*nSNPs, nSNPs, pr.i, lower.tail=F)]
  win.out[,perm_type:=ifelse(perm==0, "real","permuted")]
  win.out[,invName:=case_when(
        chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
        chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
        chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
        chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
        chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
        chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
        TRUE ~ "noInv")]

### some basic plots
  ggplot(data=win.out[nSNPs>100], aes(y=nSNPs, x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(~term)

  ggplot(data=win.out[nSNPs>25], aes(y=rnp.pr, x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(chr~term)
  ggplot(data=win.out[nSNPs>100], aes(y=wZa, x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(chr~term)
  ggplot(data=win.out[nSNPs>100], aes(y=-log10(wZa.p), x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(chr~term)

  ggplot(data=win.out[nSNPs>25], aes(y=-log10(rnp.p), x=pos_mean, group=perm, color=perm_type)) + geom_line() + facet_grid(term~chr, scales="free_y")
  ggplot(data=win.out[nSNPs>25], aes(y=rnp.pr, x=pos_mean, group=perm, color=perm_type)) + geom_line() + facet_grid(term~chr, scales="free_y")

### basic localization
  win.out.pa <- win.out[,list(wZa.pa=p.adjust(wZa.p), wZa.rnp=rank(wZa.p)/(1+length(wZa.p)), win), list(perm, term)]
  setkey(win.out, perm, term, win)
  setkey(win.out.pa, perm, term, win)

  win.out2 <- merge(win.out, win.out.pa)

  win.out2.ag <- win.out2[,list(mu=mean(wZa.pa<.05)), list(chr, invName, perm, term)]
  ggplot(data=win.out2.ag, aes(x=invName, y=mu, color=as.factor(perm==0))) + geom_point() + facet_grid(term~chr, scales="free_x")
  ggplot(data=win.out[nSNPs>25], aes(y=-log10(wZa.p), x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(chr~term)
  ggplot(data=win.out[nSNPs>100], aes(y=-log10(wZa.p), x=perm_type, group=perm_type, fill=perm_type)) + geom_boxplot() + facet_grid(chr~term)




  win.out.ag <- win.out2[perm!=0,list(minZ=min(wZa)), list(win, term)]
  setkey(win.out2, win, term)
  setkey(win.out.ag, win, term)

  win.out3 <- merge(win.out2[perm==0], win.out.ag)

  ww <- dcast(win.out3[nSNPs>25], win + chr + pos_mean + invName ~ term, value.var=c("wZa", "wZa.p", "wZa.pa", "minZ", "rnp.pr", "rnp.p"))

  ggplot(data=ww, aes(x=rnp.pr_p_lrt.x, y=rnp.pr_p_lrt.y)) + geom_point() + facet_grid(invName~chr)
  tab <- table(ww$wZa.pa_p_lrt.x < .0005 ,
               ww$wZa.pa_p_lrt.y < .0005 , ww$chr);





 ww <- dcast(win.out2[nSNPs>25], win + chr + pos_mean + invName + perm~ term, value.var=c("wZa", "wZa.p", "wZa.pa", "rnp.pr", "rnp.p"))
 ww.ag <- ww[perm!=0,list(thr.rnp.pr_p_lrt.x=quantile(rnp.pr_p_lrt.x, .9),
                          thr.rnp.pr_p_lrt.x=quantile(rnp.pr_p_lrt.y, .9)),
                      list(perm)]



 ggplot(data=ww[invName=="2Lt"], aes(x=rnp.pr_p_lrt.x, y=rnp.pr_p_lrt.y)) +
          geom_point() + facet_grid(~perm)

 ggplot(data=ww, aes(x=rnp.pr_p_lrt.x, y=rnp.pr_p_lrt.y, color=as.factor(rnp.pr_p_lrt.x>.085 & rnp.pr_p_lrt.y >.085))) +
          geom_point() + facet_grid(chr+as.factor(invName!="none")~perm) + geom_hline(yintercept=.085) + geom_vline(xintercept=.085)


  ww.ag <- ww[,list(joint=mean(rnp.pr_p_lrt.x>.085 & rnp.pr_p_lrt.y >.085),
                    pr.x=mean(rnp.pr_p_lrt.x>.085),
                    pr.y=mean(rnp.pr_p_lrt.y>.085)),
                list(perm, chr, inv=as.factor(invName!="none"))]
  ww.ag[,exp:=pr.x*pr.y]
  ww.ag[,en:=joint/exp]
  ww.ag[order(en)]





 ggplot(data=ww[invName=="2Lt"], aes(x=-log10(rnp.p_p_lrt.x), y=-log10(rnp.p_p_lrt.y))) + geom_point() + facet_grid(~perm)



  fisher.test(tab[,,1])










  ww <- dcast(win.out2[nSNs>250], win + chr + pos_mean ~ term, value.var=c("wZa", "wZa.p", "wZa.pa", "minZ"))

  tab <- table(ww$wZa.pa_p_lrt.x < .0005 & ww$wZa_p_lrt.x < ww$minZ_p_lrt.x,
               ww$wZa.pa_p_lrt.y < .0005 & ww$wZa_p_lrt.y < ww$minZ_p_lrt.y); fisher.test(tab)

  ggplot(data=ww, aes(x=-log10(wZa.pa_p_lrt.x), y=-log10(wZa.pa_p_lrt.y))) + geom_point()


  wins <- unique(ww[wZa.pa_p_lrt.x<.05 & wZa.pa_p_lrt.y <.05 & perm==0]$win)



  ggplot(data=win.out3[nSNPs>250], aes(y=-log10(wZa.p), x=pos_mean, group=perm, color=perm_type)) +
  geom_line() + facet_grid(term~chr, scales="free_y") +
  geom_point(data=win.out[perm==0][win%in%wins], aes(y=-log10(wZa.p), x=pos_mean))





#### window enrichment test
  win.out.pa <- win.out[,list(wZa.pa=p.adjust(wZa.p), wZa.rnp=rank(wZa.p)/(1+length(wZa.p)), win), list(perm, term)]
  setkey(win.out, perm, term, win)
  setkey(win.out.pa, perm, term, win)

  win.out2 <- merge(win.out, win.out.pa)

  ww <- dcast(win.out2[nSNPs>250], win + chr + pos_mean + perm ~ term, value.var=c("wZa", "wZa.p", "wZa.pa", "wZa.rnp"))


  tab <- table(ww$wZa.pa_p_lrt.x < .000005,
               ww$wZa.pa_p_lrt.y < .000005, ww$perm)


  tab <- table(ww$wZa.rnp_p_lrt.x < .05,
              ww$wZa.rnp_p_lrt.y < .05, ww$perm)

  fisher.test(tab[,,0])
  wins <- unique(ww[wZa.pa_p_lrt.x<.05 & wZa.pa_p_lrt.y <.05 & perm==0]$win)



  ggplot(data=win.out[nSNPs>250], aes(y=-log10(wZa.p), x=pos_mean, group=perm, color=perm_type)) +
  geom_line() + facet_grid(term~chr, scales="free_y") +
  geom_point(data=win.out[perm==0][win%in%wins], aes(y=-log10(wZa.p), x=pos_mean))











                  rnp.binom.p=c(binom.test(sum(seas.rnp<=pr.i),
                                           length(seas.rnp), pr.i)$p.value),
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

  save(win.out, file="~/NoCore20_NoProblems_Steep_Neg_seas_yearPop_Ran.13June2023.windows.Rdata")


  table(win.out$wZa.p<1e-10, win.out$perm)




### plot
  system("scp aob2x@rivanna.hpc.virginia.edu:~/NoCore20_NoProblems_Steep_Neg_seas_yearPop_Ran.13June2023.windows.Rdata ~/.")

  library(ggplot2)
  library(data.table)
  library(patchwork)

  load("~/NoCore20_NoProblems_Steep_Neg_seas_yearPop_Ran.13June2023.windows.Rdata")

  wZa.plot <- ggplot(data=win.out[perm!=0]) +
  geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
    geom_line(aes(x=pos_mean, y=-log10(wZa.p), group=perm), alpha=.5) +
    geom_line(data=win.out[perm==0],aes(x=pos_mean, y=-log10(wZa.p)), color="red") +
    facet_grid(~chr)+
    facet_grid(~chr) + geom_hline(aes(yintercept=-log10(min(win.out[perm!=0]$wZa.p))))


  rnp.plot <- ggplot(data=win.out[perm!=0]) +
  geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
    geom_line(aes(x=pos_mean, y=-log10(rnp.binom.p), group=perm), alpha=.5) +
    geom_line(data=win.out[perm==0],aes(x=pos_mean, y=-log10(rnp.binom.p)), color="red") +
    facet_grid(~chr) + geom_hline(aes(yintercept=-log10(min(win.out[perm!=0]$rnp.binom.p))))


window.mega <- wZa.plot / rnp.plot
ggsave(window.mega, file="~/window_mega.13June2023.pdf", h=8, w=11)



  win.out.pa <- win.out[,list(wZa.pa=p.adjust(wZa.p, "bonferroni"),
                              rnp.binom.pa=p.adjust(rnp.binom.p, "bonferroni"), win), list(perm)]




  setkey(win.out, perm, win)
  setkey(win.out.pa, perm, win)
  win.out <- merge(win.out, win.out.pa)




    wZa.plot <- ggplot(data=win.out[perm!=0]) +
    geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
      geom_line(aes(x=pos_mean, y=-log10(wZa.pa), group=perm), alpha=.5) +
      geom_line(data=win.out[perm==0],aes(x=pos_mean, y=-log10(wZa.pa), group=invName, color=invName)) +
      facet_grid(~chr)+
      facet_grid(~chr) + geom_hline(aes(yintercept=-log10(min(win.out[perm!=0]$wZa.pa))))


    rnp.plot <- ggplot(data=win.out[perm!=0]) +
    geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
      geom_line(aes(x=pos_mean, y=-log10(rnp.binom.pa), group=perm), alpha=.5) +
      geom_line(data=win.out[perm==0],aes(x=pos_mean, y=-log10(rnp.binom.pa)), color="red") +
      facet_grid(~chr) + geom_hline(aes(yintercept=-log10(min(win.out[perm!=0]$rnp.binom.pa))))


  window.mega <- wZa.plot / rnp.plot
  ggsave(window.mega, file="~/window_mega.pa.13June2023.pdf", h=8, w=11)




ggplot(data=win.out[perm==0][nSNPs>50], aes(x=-log10(wZa.p), y=wZa)) + geom_point()
ggplot(data=win.out[perm==0][nSNPs>50], aes(x=-log10(wZa.p), y=-log10(rnp.binom.p))) + geom_point()
ggplot(data=win.out[perm==0][nSNPs>50], aes(x=wZa, y=rnp.pr)) + geom_point()

win.out.ag <- win.out[nSNPs>50][rnp.pr!=Inf][,list(
                            rnp.pr=rnp.pr[perm_type=="real"],
                            rnp.perm=mean(rnp.pr[perm!=0], na.rm=T),
                            rnp.binom.p=rnp.binom.p[perm==0],
                            rnp.binom.perm=min(rnp.binom.p[perm!=0], na.rm=T)),
                        list(win, chr, pos_mean, pos_min, pos_max)]






#### old



mat <- matrix(c(681104, 13191,
                16344,    296), byrow=T, nrow=2, ncol=2)

  table(mm$p_lrt.x<1e-3, mm$perm)
  table(mm$p_lrt.y<1e-3, mm$perm)

  mm[,combo_p:=1-pchisq(-2*(log(p_lrt.x) + log(p_lrt.y)), 4)]
  table(mm$combo_p<1e-3, mm$perm, mm$chr)




  table(is.na(m$p_lrt.x), is.na(m$p_lrt.y))
  table(is.na(m$p_lrt.x), m$p_lrt.y<1e-3)
  table(m[chr=="2R"]$p_lrt.x < 1e-3,  m[chr=="2R"]$p_lrt.y <1e-3, m[chr=="2R"]$perm)
  table(m$p_lrt.x < 1e-3, m$perm,  m$chr)
  table(m$p_lrt.y < 1e-3,  m$perm, m$chr)

  m.ag <- m[,list(or=fisher.test(p_lrt.x<.05, p_lrt.y<.05)$estimate), list(perm, chr)]
  m.ag[order(or)]

  m[,.N,perm]
