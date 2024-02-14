# ijob -A biol4559-aob2x -c10 -p largemem --mem=240G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(tidyr)
  library(dplyr)

### load data
  load("~/dest2_glm_baypass_annotation.Rdata") ### made by `DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/collectBayPass.R`

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
  m[,XtXst.rnp:=rank(-XtXst_median)/(length(XtXst_median)+1)]
  m[XtXst.rnp<1e-5]

  m[,C2.rnp:=rank(-C2_std_median)/(length(C2_std_median)+1)]
  m[C2.rnp<1e-5]


  m[,C2.rnp:=rank(-C2_std_median)/(length(C2_std_median)+1)]
  m[C2.rnp<1e-5]

  m[,prod:=-log(XtXst.rnp*C2.rnp)]

  perm <- foreach(i=0:100, .combine="rbind")%do%{
    if(i==0) return(data.table(prod=-log10(m$XtXst.rnp * (m$C2.rnp)), i=i))
    if(i!=0) return(data.table(prod=-log10(m$XtXst.rnp * sample(m$C2.rnp)), i=i))
  }

  ggplot() +
  geom_density(data=perm[i>0 & i<10], aes(prod, group=i), alpha=.5) +
  geom_density(data=perm[i==0], aes(prod, group=i), color="red") +
  geom_density(data=data.table(prod=rchisq(10000, 2)), color="red")


  mean(quantile(-log10(m$XtXst.rnp * (m$C2.rnp)), .999999) > thr)
  mean(quantile(-log10(m$XtXst.rnp * (m$C2.rnp)), .999999) / thr)

  ggplot(data=m[prod>5], aes(x=pos, y=prod)) + geom_point(size=.5, alpha=.5) + facet_grid(~chr)

  ggplot(data=m, aes(prod)) + geom_density() +
  geom_density(data=data.table(prod=rchisq(10000, 2)), color="red")
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
    win.tmp[,C2.Z:=qnorm(cont.p, 0, 1)]
    win.tmp[,xtx.Z:=qnorm(xtx.p, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]

    win.tmp[,het:=2*af*(1-af)]



    tab1 <- table(win.tmp$XtXst.rnp<.05, win.tmp$rnp<.05)
    tab2 <- table(win.tmp$XtXst.rnp<.05, win.tmp$C2.rnp<.05)
    tab3 <- table(win.tmp$rnp<.05, win.tmp$C2.rnp<.05)

    fet1 <- fisher.test(tab1)
    fet2 <- fisher.test(tab2)
    fet3 <- fisher.test(tab3)

    pr.i =c(0.05)

      win.out <- win.tmp[,list(chr=chr[1], pos_mean = mean(pos), pos_min = min(pos), pos_max = max(pos),
                    win=win.i, pr=pr.i,
                    glm.rnp.pr=c(mean(rnp<=pr.i)),
                    glm.rnp.binom.p=c(binom.test(sum(rnp<=pr.i),
                                             length(rnp), pr.i, alternative="greater")$p.value),
                   xtx.rnp.pr=c(mean(XtXst.rnp<=pr.i)),
                   xtx.rnp.binom.p=c(binom.test(sum(XtXst.rnp<=pr.i),
                                            length(XtXst.rnp), pr.i, alternative="greater")$p.value),
                    c2.rnp.pr=c(mean(C2.rnp<=pr.i)),
                    c2.rnp.binom.p=c(binom.test(sum(C2.rnp<=pr.i),
                                             length(C2.rnp), pr.i, alternative="greater")$p.value),
                    wZa=sum(het*Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                    wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                    C2.wZa=sum(het*C2.Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                    C2.wZa.p=pnorm(sum(het*C2.Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                    xtx.wZa=sum(het*xtx.Z, na.rm=T)/(sqrt(sum(het^2, na.rm=T))),
                    xtx.wZa.p=pnorm(sum(het*xtx.Z)/(sqrt(sum(het^2))), lower.tail=T, log.p=T),
                    rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                    rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                    min.p.lrt=min(p_lrt),
                    min.rnp=min(rnp),
                    nSNPs = .N,
                    sum.rnp=sum(rnp<=pr.i),
                    glm.genes=paste(unique(gene[rnp<=pr.i]), collapse=";"),
                    xtx.genes=paste(unique(gene[XtXst.rnp<=pr.i]), collapse=";"),
                    c2.genes= paste(unique(gene[C2.rnp<=pr.i]), collapse=";"),
                    combo.genes=paste(unique(gene[rnp<=pr.i & C2.rnp<=pr.i & XtXst.rnp<=pr.i]), collapse=";"),
                    delux.combo.genes=paste(unique(gene[rnp<=pr.i & C2.rnp<=pr.i & XtXst.rnp<=pr.i & col=="missense_variant"]), collapse=";"),
                    xtx_glm.or=fet1$estimate, xtx_c2.or=fet2$estimate, glm_c2.or=fet3$estimate,
                    xtx_glm.p=fet1$p.value, xtx_c2.p=fet2$p.value, glm_c2.p=fet3$p.value,
                    xtx_median=median(XtXst_median, na.rm=T), c2_median=median(C2_std_median, na.rm=T),
                    xtx_95=quantile(XtXst_median, .95, na.rm=T), c2_95=quantile(C2_std_median, .95, na.rm=T)),]

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

  save(win.out, file="~/XtX_C2_glm.windows.Rdata")


  table(win.out$wZa.p<1e-10, win.out$perm)





  cSNPS.ag <- cSNPs[,list(.N), list(col, gene)]

  cSNPS.ag[col=="missense_variant"][order(N)]











### overlap plots
  xtx_c2_overlap.plot <-
    ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
    geom_line(data=win.out,aes(x=pos_mean, y=log2( xtx_c2.or)), color="red") +
    facet_grid(~chr, scales="free_x") + ylab("log2(xtx_c2 odds)")


### nSNps and spread
  ggplot(data=win.out, aes(x=nSNPs, y=c2_median)) + geom_point() + geom_vline(xintercept=100)


table(win.out$)


win.out[C2.wZa.p< -5][chr=="3L"][pos_mean>10e6 & pos_mean<12e6]
win.out[xtx.wZa.p< -500][C2.wZa.p < -12][wZa.p< -75]

ggplot(data=win.out, aes(x=log2(xtx_c2.or), y=-log10(xtx_c2.p))) + geom_point()
win.out[which.min(xtx_c2.p)]
  xtx_overlap <- ggplot(data=win.out) +
   geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
     geom_line(data=win.out,aes(x=pos_mean, y=-log10(xtx_c2.p)*sign(log2(xtx_c2.or))), color="red") +
     facet_grid(~chr) + ylab("XtX-c2 overlap p-value")

   xtx_median.plot <- ggplot(data=win.out) +
   geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
     geom_line(data=win.out,aes(x=pos_mean, y=xtx_median), color="red") +
     facet_grid(~chr)


layout <- "
A
B
C"

c2.wza.plot + xtx_median.plot + xtx_overlap + plot_layout(design=layout)



win.out[-log10(C2.wZa.p)>15][1]


  glm.rnp.plot <- ggplot(data=win.out) +
  geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
    geom_line(data=win.out,aes(x=pos_mean, y=-log10(glm.rnp.binom.p), group=pr, color=as.factor(pr))) +
    facet_grid(pr~chr)

  xtx.rnp.plot <- ggplot(data=win.out) +
  geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
    geom_line(data=win.out,aes(x=pos_mean, y=-log10(xtx.rnp.binom.p)), color="red") +
    facet_grid(pr~chr)

  c2.rnp.plot <- ggplot(data=win.out) +
  geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
    geom_line(data=win.out,aes(x=pos_mean, y=-log10(C2.wZa.p)), color="red") +
    facet_grid(~chr)

    xtx_overlap <- ggplot(data=win.out) +
   geom_vline(data=data.table(chr=c("3R", "3L"), pos_mean=c(13223000, 3317423)), aes(xintercept=pos_mean), color="green") +
     geom_line(data=win.out,aes(x=pos_mean, y=-log10(xtx_c2.p)*sign(log2(xtx_c2.or))), color="red") +
     facet_grid(~chr)


   c2.rnp.plot <- ggplot(data=win.out) +
   geom_vline(data=data.table(chr=c("3R", "3L"), pos_mean=c(13223000, 3317423)), aes(xintercept=pos_mean), color="green") +
     geom_line(data=win.out,aes(x=pos_mean, y=-log10(c2.rnp.binom.p), group=pr, color=as.factor(pr))) +
     facet_grid(pr~chr)


   wza.plot <- ggplot(data=win.out) +
   geom_vline(data=data.table(chr=c("3R", "3L"), pos_mean=c(13223000, 3317423)), aes(xintercept=pos_mean), color="green") +
   geom_line(data=win.out,aes(x=pos_mean, y=-log10(wZa.p), group=pr, color=as.factor(pr))) +
   facet_grid(pr~chr)


  xtx_median / glm.rnp.plot








    ggplot(data=win.out, aes(x=xtx_median, y=-log10(glm.rnp.binom.p), color=xtx_glm.or)) + geom_point()

    win.out[which.max(-log10(xtx_c2.p)*sign(log2(xtx_c2.or)))]
    win.out[-log10(wZa.p)>20]
     ggplot(data=win.out[nSNPs>200]) +
    geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
      geom_line(data=win.out,aes(x=pos_mean, y=c2_95), color="red") +
      facet_grid(~chr)


win.out[chr=="3R"][which.min(glm.rnp.binom.p)]


window.mega <- wZa.plot / rnp.plot
ggsave(window.mega, file="~/window_mega.13June2023.pdf", h=8, w=11)


win.out[xtx.rnp.binom.p]








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

win.out.real <- win.out[perm==0]
win.out.perm <- win.out[perm!=0][,list(rnp.binom.p.perm=min(rnp.binom.p), rnp.pr.perm=median(rnp.pr),
                                      wZa.p.perm=min(wZa.p), wZa.perm=mean(wZa)), win]
win.out.ag <- merge(win.out.real, win.out.perm, by="win", all=T)
win.out.ag[,rnp.en:=rnp.pr/rnp.pr.perm]
win.out.ag[,wZa.d:=wZa - wZa.perm]

ggplot(data=win.out.ag) +
geom_vline(data=data.table(chr="3R", pos_mean=13223000), aes(xintercept=pos_mean), color="green") +
  geom_line(data=win.out.ag[perm==0],aes(x=pos_mean, y=wZa.d)) +
  facet_grid(~chr)











win.out[win==2163]
paste(
paste(
  win.out[wZa.pa<.05 & wZa.p<=min(win.out[perm!=0]$wZa.p)]$chr, ":",
  win.out[wZa.pa<.05 & wZa.p<=min(win.out[perm!=0]$wZa.p)]$pos_min, "-",
  win.out[wZa.pa<.05 & wZa.p<=min(win.out[perm!=0]$wZa.p)]$pos_max, sep=""),
collapse=";")


2L:18705770-18805560
2L:19056063-19155754
2L:22507054-22604751
2L:22557817-22646654
2R:8426164-8525994
2R:10026970-10125236
2R:13476159-13575942
2R:14776341-14875606
2R:19776125-19875871
3L:1921175-2020903
3L:11720962-11820923
3L:19571411-19670889
3L:19621269-19708428
3L:20170947-20270562
3L:20221019-20320938
3L:20271499-20370642
3L:21821069-21920452
3L:21870985-21970634
3R:5836746-5936414
3R:6286865-6384777
3R:10086544-10186072
3R:10136542-10231655
3R:10186945-10286340
3R:10237670-10336250
3R:12536489-12636185
3R:21086486-21186446
3R:21136679-21236359
3R:29036833-29136414
2R:5526482-5621991







2L:1005866-1104670
2L:1056969-1155526
2L:2606211-2705721
2L:6305815-6405723
2L:6505812-6605566
2L:9155769-9255533
2L:9905795-10005703
2L:12356071-12455524
2L:12406026-12505664
2L:13356119-13449465
2L:18656323-18755752
2L:18705770-18805560
2L:19005784-19105532
2L:19056063-19155754
2L:19106171-19205677
2L:19506514-19605719
2L:19555798-19652528
2L:20555874-20653149
2L:22507054-22604751
2L:22557817-22646654
2R:8426164-8525994
2R:8976462-9075956
2R:9326502-9426113
2R:10026970-10125236
2R:13076217-13175995
2R:13476159-13575942
2R:14776341-14875606
2R:14926147-15026121
2R:15626610-15726043
2R:16226136-16326109
2R:16276212-16376080
2R:17476237-17576082
2R:18376310-18475992
2R:19727323-19826091
2R:19776125-19875871
2R:20876287-20975942
2R:21026185-21126108
2R:21076169-21176104
2R:21126246-21226080
2R:21176166-21276122
2R:21226210-21326052
2R:21326802-21425888
2R:21626294-21725954
2R:21676266-21776120
2R:21876142-21976108
2R:22176287-22276048
2R:24576185-24675405
3L:1771069-1870724
3L:1921175-2020903
3L:3271121-3370713
3L:3321944-3420595
3L:3871080-3970537
3L:4221459-4320172
3L:4270947-4370798
3L:4771131-4870727
3L:5124630-5220903
3L:9571011-9670648
3L:9621237-9720932
3L:9671053-9770911
3L:10120992-10220822
3L:10170951-10270804
3L:10220983-10320899
3L:11720962-11820923
3L:18021058-18120931
3L:18523164-18620788
3L:19421215-19520854
3L:19471646-19570654
3L:19521150-19620632
3L:19571411-19670889
3L:19621269-19708428
3L:20170947-20270562
3L:20221019-20320938
3L:20271499-20370642
3R:5686580-5785403
3R:5836746-5936414
3R:6736683-6836180
3R:6786483-6886430
3R:10086544-10186072
3R:10136542-10231655
3R:10237670-10336250
3R:10686660-10785369
3R:10736658-10836369
3R:11286475-11379063
3R:11336720-11436320
3R:12536489-12636185
3R:13186559-13286443
3R:14686466-14786423
3R:18286467-18386429
3R:18336512-18435636
3R:18736472-18836394
3R:21086486-21186446
3R:21136679-21236359
3R:22236626-22336430
3R:23936638-24036209
3R:25286505-25386399
3R:28686509-28785750
3R:28736760-28836343
3R:29036833-29136414
3R:29586496-29686421
3R:29636584-29736437
2R:12176201-12276087



  win.out.ag <- win.out[,list()]
