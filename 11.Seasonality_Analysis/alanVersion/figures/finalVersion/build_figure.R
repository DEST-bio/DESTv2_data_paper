### library
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load data
  setwd("~/dest2_seasonality")

  ### A, B, C: this file loads the output of the WZA on the un-permuted data. WZA for XtXst pvalue (Chisq based), Contrast pvalue, and GLM pvalue
    ### downloaded
    #system("scp aob2x@rivanna.hpc.virginia.edu:~/XtX_C2_glm.windows.pod_v2.Rdata ~/dest2_seasonality/XtX_C2_glm.windows.pod_v2.Rdata")
    load(file="XtX_C2_glm.windows.pod_v2.Rdata") ### `win.out`

  ### B:
    load(file="destv2_seasonality_C2.perm.Rdata") ### `win.all.ag2`

  ### C: this file contains the GLM permutations used to draw significance threshold for panel C.
    ### downloaded
    load(file="destv2_seasonality_perm.Rdata") ### `win.all.ag2`
    thrs <- 1-c(.95, .99, .995, .999)
    glm.permuation.wZA.thresholds <- win.all[perm!=0,list(value=quantile(wZa.p, thrs), thr=thrs), list(perm)]
    glm.perm.ag <- glm.permuation.wZA.thresholds[,list(value=min(value)), list(thr)]
    #win.all.ag2[,pos_mean:=pos_min/2 + pos_max/2]

  ### D: this file contains the seasonal pairs. reference this file: `DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures/mapFigure.R`
    samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

    seasonal.sets <- samps[Recommendation=="Pass"][continent%in%c("Europe", "Asia")][set!="DrosRTEC"][,list(.N,
                                                                   season=c("spring", "fall"),
                                                                   sampleId=c(sampleId[which.min(jday)], sampleId[which.max(jday)]),
                                                                   jday=c(jday[which.min(jday)], jday[which.max(jday)]),
                                                                   diff=c(jday[which.min(jday)]-jday[which.max(jday)]),
                                                                   lat=lat[which.min(jday)]),
                                                                list(country, city, locality, year, loc_rep, cluster2.0_k5)][N>1][abs(diff)>30]


  ### E: this file should contain what we want
    load(file="enrichment.NoCore20_seas_europe_yearPop_Ran.Rdata")

  ### F: this file contains the detailed SNP level metrics that allow us to build the hexbin plot. Also contains the thresholds for the POD analysis of C2 and XTX
    ### downloaded
    #system("scp aob2x@rivanna.hpc.virginia.edu:~/dest2_glm_baypass_annotation_pod.podOutputToo_v2.Rdata ~/dest2_seasonality/.")
    load(file="dest2_glm_baypass_annotation_pod.podOutputToo_v2.Rdata") ###thrs.ag contains C2 thresholds
    load("glm.perm.thr.ag.Rdata")

### build panels
  nSNP_thr <- 150

  ### XTX wza
    xtx.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-xtx.wZa.p), color="black") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(xtx wZa p)") + xlab("Pos (Mb)") + theme_bw()

    xtx.wza.pod.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-xtx.wZa.pod.p), color="black") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(xtx wZa p)") + xlab("Pos (Mb)") + theme_bw()

    xtx.wza.plot /xtx.wza.pod.plot

    ggplot(data=win.out, aes(x=xtx.wZa.pod, y=xtx.wZa)) + geom_point()
    win.out[chr=="2L"][which.min(xtx.wZa.pod.p)]

  ### C2 WZA
    C2.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-C2.wZa.p), color="black") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p)") + xlab("Pos (Mb)") + theme_bw()

    C2.wza.pod.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-C2.wZa.pod.p), color="black") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p pod)") + xlab("Pos (Mb)") + theme_bw() +
    geom_hline(yintercept=-log10(.0001/dim(win.out)[1]))

    C2.wza.plot /C2.wza.pod.plot
    ggplot(data=win.out, aes(x=C2.wZa.pod, y=C2.wZa)) + geom_point()

  ### GLM WzA
    GLM.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-wZa.p), color="black") +
    geom_hline(data=glm.perm.ag[4], aes(yintercept=-value, color=as.factor(thr))) +
    facet_grid(~chr, scales="free_x") + ylab("-log10(GLMER wZa p)") + xlab("Pos (Mb)") + theme_bw() + theme(legend.position="none")

  ### Collection Plot
    ss.ag <- seasonal.sets[,list(lat=mean(lat)), list(locality, city, country)]
    ss.ag <- ss.ag[order(lat)]
    ss.ag[,ypos:= seq(from=min(lat), to=max(lat), length.out=length(lat))]
    ss.ag

    collection_fig <- ggplot(data=seasonal.sets, aes(x=jday, y=lat, group=interaction(locality, year))) +
    geom_hline(data=seasonal.sets[,list(lat=mean(lat)), list(locality)], aes(yintercept=lat), color="grey", linetype="dashed") +
    geom_line() +
    geom_point() + ylab("Latitude") + xlab("Julian Day") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~year)

    line_fig <- ggplot(data=ss.ag) +
    geom_segment(aes(x=0, xend=1, y=lat, yend=ypos), color="grey", linetype="dashed") +
    theme_void() + xlim(0,1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


    label_fig <- ggplot(data=ss.ag) +
    geom_text(data=ss.ag, aes(label=paste(city, ",\n", country, sep=""), x=.45, y=ypos), size=2.5, hjust = 0) +
    theme_void() + xlim(0,1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ### GLM enrichment relative to permutations
    oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, min_p, max_p, model_features, pops)]

    oo.ag.ag <- oo.ag[,list(N=N[perm==0],
                      N.perm.mean=mean(N[perm!=0]),
                      N.perm.lci=quantile(N[perm!=0], .005),
                      N.perm.uci=quantile(N[perm!=0], .995)),
                  list(chr, min_p, max_p, pops, model_features)]


    glm.real_vs_perm.plot <- ggplot(data=oo.ag.ag) +
    geom_ribbon(aes(x=-log10(max_p), ymin=(N.perm.uci), ymax=N.perm.lci), alpha=.5, fill="grey") +
    geom_line(aes(x=-log10(max_p), y=(N.perm.mean)), alpha=.5) +
    geom_line(aes(x=-log10(max_p), y=(N)), color="red") +
    xlab("-log10(p_lrt)") + theme_bw()

  ### SNP level intersection between GLM and C2
    fisher.test(table(-log10(m$p_lrt)>3.5, -log10(m$cont.pod.p)>3.5))
    fisher.test(table(-log10(m$p_lrt)>3.5, (m$XtXst_median)>250))

    lab <- paste("Odds Ratio=", round(as.numeric(fisher.test(table(-log10(m$p_lrt)>-log10(glm.perm.thr.ag[thr==1-.999]$glm_thr),
                                                                          m$C2_std_median>thrs.ag[thr==.999]$C2_thr))$estimate)), sep="")


    glm_c2 <- ggplot(m) +
    geom_hex(aes(x=-log10(p_lrt), y=C2_std_median)) +
    geom_vline(xintercept=thrs.ag[thr==.999]$C2_thr, color="red") +
    geom_hline(yintercept=-log10(glm.perm.thr.ag[thr==1-.999]$glm_thr), color="red") +
    #geom_text(aes(x=6, y=20, label=lab)) +
    theme_bw() + theme(legend.position="none")



### combined plot

  layout <- "
  AAAAAABBCCCCCCCCCCCC
  AAAAAABBDDDDDDDDDDDD
  AAAAAABBEEEEEEEEEEEE
  AAAAAABBFFFFFFGGGGGG
  AAAAAABBFFFFFFGGGGGG"


  layout <- "
  AAAAAABBCCCCCCCCCCCC
  AAAAAABBCCCCCCCCCCCC
  AAAAAABBDDDDDDDDDDDD
  AAAAAABBDDDDDDDDDDDD
  AAAAAABBEEEEEEEEEEEE
  AAAAAABBEEEEEEEEEEEE
  FFFFFFFFFFFGGGGGGGGG
  FFFFFFFFFFFGGGGGGGGG
  "

  mega <- collection_fig + label_fig + xtx.wza.plot + C2.wza.plot + GLM.wza.plot + glm.real_vs_perm.plot + glm_c2 + plot_layout(design=layout)

### version 2
  layout <- "
  AAAAABHHH
  AAAAABHHH
  AAAAABHHH
  AAAAABHHH
  AAAAABHHH
  CCCDDDHHH
  CCCDDDHHH
  CCCDDDHHH
  EEEEEEEEE
  EEEEEEEEE
  FFFFFFFFF
  FFFFFFFFF
  GGGGGGGGG
  GGGGGGGGG"
  mega <- collection_fig + line_fig + glm.real_vs_perm.plot + glm_c2 + xtx.wza.pod.plot + C2.wza.pod.plot + GLM.wza.plot + label_fig + plot_layout(design=layout)


  ggsave(mega, file="seasonal_mega.pod.pdf", w=7, h=12)





### supplemental
  dim(m)
  dim(xtx.ag)

  ### XtX observed vs POD
    m[,XtXst_rank:=NA]
    xtx.ag[,XtXst_rank:=NA]

    m[!is.na(XtXst_median),XtXst_rank:=rank(XtXst_median)]
    xtx.ag[!is.na(XtXst_median),XtXst_rank:=rank(XtXst_median)]
    xtx.ag[which.max(XtXst_rank)]

    xtx.dist <- merge(m[,c("XtXst_rank", "XtXst_median", "xtx_neglogp_median"), with=F], xtx.ag[,c("XtXst_rank", "XtXst_median", "xtx_neglogp_median"), with=F], by="XtXst_rank")

    qq_xtx <- ggplot(data=xtx.dist, aes(y=XtXst_median.x, x=XtXst_median.y)) + geom_point()
    ggsave(qq_xtx, file="~/qq_xtx_v2.png")

    ggplot(data=m, aes(XtXst_median), fill="red") + geom_density() +
    geom_density(data=xtx.ag, aes(XtXst_median), alpha=.75, fill="black")

    ggplot(data=m, aes(C2_std_median), color="red") + geom_density() +
    geom_density(data=cont.ag, aes(C2_std_median), alpha=.75, fill="black")

    fisher.test(rbind(table(m$C2_std_median>25), table(cont.ag$C2_std_median>25)))


    ggplot(data=xtx.ag, aes(x=XtXst_median, y=xtx_neglogp_median)) + geom_hex()
    ggplot(data=m, aes(x=XtXst_median, y=xtx_neglogp_median)) + geom_hex()
    ggplot(data=m, aes(x=xtx_neglogp_median, y=-log10(XtXst.pod.p))) + geom_hex()
    ggplot(data=m, aes(XtXst.pod.p)) + geom_histogram()


    fisher.test(table(m$chr=="2L" & m$invName.x!="noInv", m$XtXst_median<75))


    glm_c2 <- ggplot(m) +
    geom_hex(aes(x=-log10(cont.pod.p), y=C2_std_median)) +
    #geom_text(aes(x=6, y=20, label=lab)) +
    theme_bw() + theme(legend.position="none")


  t(table(m$XtXst.pod.p<.5e-5, m$col))
  (559/34295)/(1673/316573)
  m[col=="stop_gained"][XtXst.pod.p<5e-5]
  m[col=="stop_gained"][XtXst.pod.p<5e-5]


  t(table(m$cont.pod.p<.5e-5, m$col))
  m[grepl("stop_gained", col)][cont.pod.p<5e-5]

(145/34709)/(454/115704)


fisher.test(table(m$cont.pod.p<.5e-5, m$XtXst.pod.p<.5e-5))
fisher.test(table(m$cont.pod.p<.5e-5 & m$XtXst.pod.p<.5e-5, m$col=="missense_variant"))
fisher.test(table(m$cont.pod.p<.5e-5, m$col=="missense_variant"))
fisher.test(table(m$XtXst.pod.p<.5e-5, m$col=="missense_variant"))

m[cont.pod.p<.5e-5 & XtXst.pod.p<.5e-5 & col=="missense_variant"]

### C2
ggplot(data=m, aes(C2_std_median)) + geom_density() +
geom_density(data=cont.ag, aes(C2_std_median), fill="red")
