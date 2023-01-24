### libraries
  library(data.table)
  library(doMC)
  registerDoMC(4)
  library(ggplot2)
  library(patchwork)

### get PCR dup rates for all libraries
#   fns <- system("find /project/berglandlab/DEST/dest_mapped/ -name '*duplicates_report.txt'", intern=T)
#   pcr <- foreach(fn=fns)%dopar%{
#     #fn <- fns[1]
#     message(fn)
#     tmp <- fread(fn, skip="LIBRARY", nrows=1)
#     data.table(pcrDup=tmp$PERCENT_DUPLICATION, READ_PAIRS_EXAMINED=tmp$READ_PAIRS_EXAMINED, sampleId#=tstrsplit(fn, "/")[[7]])
#   }
#   pcr <- rbindlist(pcr)
#   save(pcr, file="./pcr_destV2.Rdata")

  # scp aob2x@rivanna.hpc.virginia.edu:~/pcr_destV2.Rdata ~/.
### load precomputed data

  ### pcr dup
    # load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/1.Quality_Control/reseq/dup.rate.DEST2.Rdata")
    # dup <- as.data.table(dup.rate)
    #
    #
    load("./pcr_destV2.Rdata")
    dup <- pcr

    dup[sampleId=="AU_Que_Inn_1_NA", sampleId:="AU_Que_Inn_NA-MM-DD"]
    dup[sampleId=="AU_Vic_Yer_1_NA", sampleId:="AU_Vic_Yer_NA-MM-DD"]
    dup[sampleId=="GB_Fif_Dai_1_NA", sampleId:="GB_Fif_Dai_2020-MM-DD"]
    dup[sampleId=="GR_Kri_Ira_1_NA", sampleId:="GR_Kri_Ira_2019-MM-DD"]
    dup[sampleId=="RS_Kol_Sla_1_NA", sampleId:="RS_Kol_Sla_2021-MM-DD"]
    dup[sampleId=="UA_Lvi_Dro_1_2014-08-24", sampleId:="UA_Lviv_Dro_1_2014-08-24"]
    dup[sampleId=="UA_Lvi_Dro_1_2015-09-18", sampleId:="UA_Lviv_Dro_1_2015-09-18"]
    dup[sampleId=="UA_Lvi_Dro_1_2016-09-12", sampleId:="UA_Lviv_Dro_1_2016-09-12"]


  ### simulans
    load("./sim.contam.joint.Rdata")
    sim <- as.data.table(sim.cont)

  ### coverage
    load("~/depth_per_contig.Rdata")
    o[,sampleId:=gsub(".original.bam", "", samp)]
    o[,chr:=gsub("sim_", "", chrom)]
    o.ag <- o[,list(cov=mean(mean[!grepl("sim", chrom)]), sim_alan=mean(mean[grepl("sim", chrom)]/sum(mean))), list(sampleId, chr)]
    o.ag.ag <- o.ag[,list(cov=mean(cov), sim_alan=mean(sim_alan)), list(sampleId)]

    o.ag.ag[sampleId=="AU_Que_Inn_1_NA", sampleId:="AU_Que_Inn_NA-MM-DD"]
    o.ag.ag[sampleId=="AU_Vic_Yer_1_NA", sampleId:="AU_Vic_Yer_NA-MM-DD"]
    o.ag.ag[sampleId=="GB_Fif_Dai_1_NA", sampleId:="GB_Fif_Dai_2020-MM-DD"]
    o.ag.ag[sampleId=="GR_Kri_Ira_1_NA", sampleId:="GR_Kri_Ira_2019-MM-DD"]
    o.ag.ag[sampleId=="RS_Kol_Sla_1_NA", sampleId:="RS_Kol_Sla_2021-MM-DD"]
    o.ag.ag[sampleId=="UA_Lvi_Dro_1_2014-08-24", sampleId:="UA_Lviv_Dro_1_2014-08-24"]
    o.ag.ag[sampleId=="UA_Lvi_Dro_1_2015-09-18", sampleId:="UA_Lviv_Dro_1_2015-09-18"]
    o.ag.ag[sampleId=="UA_Lvi_Dro_1_2016-09-12", sampleId:="UA_Lviv_Dro_1_2016-09-12"]



  ### samps
    samps <- fread("/scratch/yey2sn/DESTv2/populationInfo/dest_v2.samps_19Jan2023.csv")

  ### merge all together
    dim(samps)
    dim(dup)
    dim(pcr)

    m <- merge(dup, o.ag.ag, all=T)
    dim(m)
    m <- merge(samps[set!="dgn"], m, all=T)
    dim(m)
    table(m[is.na(sim_alan)]$set)
    table(is.na(m$sim_alan), m$low_qual)
    m[is.na(sim_alan) & low_qual==F]$sampleId
    m[is.na(locality)]$sampleId

    m[,effCov:=(2*nFlies*cov*(1-sim_alan))/(2*nFlies*(1-sim_alan)+cov)]

### basic plot

  simrate_plot <-
  ggplot(data=m, aes(x=sim_alan)) +
  geom_histogram() +
  facet_wrap(~set, nrow=1) +
  scale_x_log10() +
  annotation_logticks(sides="b", outside=T) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(vjust = -2)) +
  xlab("log10(simulans rate)") + ylab("(N)") + ggtitle("D. simulans contamination")

  pcr_plot <-
  ggplot(data=m, aes(x=pcrDup)) +
  geom_histogram() +
  facet_wrap(~set, nrow=1) +
  scale_x_log10() +
  annotation_logticks(sides="b", outside=T) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(vjust = -2), axis.title.x=element_text(vjust=-2)) +
  xlab("log10(PCR rate)") + ylab("(N)") + ggtitle("PCR dup rate")

  cov_plot <-
  ggplot(data=m, aes(x=effCov)) +
  geom_histogram() +
  facet_wrap(~set, nrow=1) +
  scale_x_log10() +
  annotation_logticks(sides="b", outside=T) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(vjust = -2), axis.title.x=element_text(vjust=-2)) +
  xlab("log10(effective coverage)")+ ylab("(N)") + ggtitle("Effective coverage")

  layout <-"
  A
  B
  C"

  seqStats.mega <-
  simrate_plot + pcr_plot + cov_plot +
  plot_annotation(tag_levels="A") + plot_layout(design=layout)
seqStats.mega
### seaonsal set
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/1.Quality_Control/reseq/frost.final.set.Rdata")
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/1.Quality_Control/reseq/DEST2.seasonals.plusCore20.flip.met.Rdata")

  seas <- as.data.table(DEST2.seasonals.plusCore20.flip.met)
  frost <- as.data.table(samps.frost.final.set)

  sf <- rbind(data.table(group="seas", sampleId=seas$sampleId), data.table(group="frost", sampleId=frost$sampleId))
  setkey(sf, sampleId)
  sf <- sf[!duplicated(sf)]

  m2 <- merge(m, sf, all=T)
  m2[is.na(locality)]
  m2[is.na(group), group:="none"]


### benefit of double coverage?
  m2[,dbl_effCov:=(2*nFlies*cov*(1-sim_alan)*2)/(2*nFlies*(1-sim_alan)*2+cov)]
  m2[,rr:=dbl_effCov/effCov]

### which are the new localities
  m2.ag <- m2[,.N,locality]

  m2 <- merge(m2, m2.ag, by="locality")
  m2[low_qual==T & N==1]


  m2[set=="DrosEU_3" & group!="none" & effCov<25 & sim_alan<.1 & effCov>15, reseq:="dbl_cov_seasonal"]
  m2[set=="DrosEU_3" & effCov<25 & sim_alan<.1 & N==1, reseq:="dbl_cov_spatial"]


  m2[set=="DrosEU_3" & group!="none" & low_qual==T & locality!="US_Vir_Cha", reseq:="need_pair_seasonal"]
  m2[set=="DrosEU_3" & N==1 & low_qual==T & locality!="US_Vir_Cha", reseq:="new_locality_spatial"]

  m2[sampleId=="MA_Tan_Lar_1_2021-09-06", reseq:="need_pair_seasonal_questionableSim"]
  write.csv(m2[!is.na(reseq), c("sampleId", "SequencingId", "reseq"), with=F], file="~/reseq.csv", quote=F)
table(m2$reseq)

  sum(table(m2$reseq))


  m2[!is.na(reseq), c("sampleId", "SequencingId", "reseq"), with=F]


  m2[locality%in%m2[reseq=="need_pair_seasonal"]$locality]$sim_alan

### output figures
  dbl.plot <- ggplot(data=m2[set=="DrosEU_3"], aes(x=effCov, y=dbl_effCov, color=reseq)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) + ylim(0, 42) + xlim(0, 42)


### plots








  ggplot(data=m2, aes(rr)) + geom_histogram()

  table(is.na(m2$group), m2$set)
  table(is.na(m2$group), m2$low_qual)







  m[sampleId%in%unique(seas$sampleId, frost$sampleId),seasonalSet]



  1: DE_Bay_Mun_1_2014-06-19  DE_Mun_14_31 DE_Bay_Mun 48.18 11.61    Europe
   2: DE_Bay_Mun_1_2014-09-03  DE_Mun_14_32 DE_Bay_Mun 48.18 11.61    Europe
   3: DE_Bay_Mun_1_2015-06-01  DE_Mun_15_45 DE_Bay_Mun 48.18 11.61    Europe
   4: DE_Bay_Mun_1_2015-09-01  DE_Mun_15_46 DE_Bay_Mun 48.18 11.61    Europe
   5: DE_Bay_Mun_1_2016-06-25  DE_Mun_16_12 DE_Bay_Mun 48.18 11.61    Europe
   6: DE_Bay_Mun_1_2016-09-01  DE_Mun_16_13 DE_Bay_Mun 48.18 11.61    Europe
   7: DE_Bay_Mun_1_2017-06-20          <NA> DE_Bay_Mun 48.18 11.61    Europe

   8: DE_Bay_Mun_1_2017-09-03          <NA> DE_Bay_Mun 48.18 11.61    Europe
   9: DE_Bay_Mun_1_2018-09-09          <NA> DE_Bay_Mun 48.18 11.61    Europe
  10: DE_Bay_Mun_1_2019-06-22          <NA> DE_Bay_Mun 48.18 11.61    Europe
  11: DE_Bay_Mun_1_2019-09-02          <NA> DE_Bay_Mun 48.18 11.61    Europe
  12: DE_Bay_Mun_1_2020-09-05          <NA> DE_Bay_Mun 48.18 11.61    Europe
  13: DE_Bay_Mun_1_2021-06-25          <NA> DE_Bay_Mun 48.18 11.61    Europe


  ggplot(data=m, aes(x=READ_PAIRS_EXAMINED*(1-sim_alan)*(1-pcrDup), y=effCov)) +
  geom_point(aes(color=sim_alan>.25)) +
  facet_wrap(~set, scale="free")




  hist(dup$PERCENT_DUPLICATION)
