### libraries
  library(data.table)
  library(patchwork)
  library(ggplot2)

### load data
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/first_batch_qc/qc.Rdata")
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/first_batch_qc/depth_per_contig.Rdata")
  setnames(o, "samp", "sampleId")

### merge
  dt <- merge(qual.dt, o, by="sampleId", allow.cartesian=T)
  dt.ag <- dt[,list(mu=mean(mean)), list(chrom)]

  dt <- merge(dt, dt.ag, by="chrom")
  dt.ag <- dt[,list(delta=(mean-mu)/mu, mean=mean), list(chrom, sampleId)]

  dt.ag[sampleId=="VA_ch_21_Aug_2020"][order(delta)][grepl("sim_2L|2L", chrom)]


### simulans rate
  dt.sim <- dt[,list(mel_rate=bases[chrom=="2L"]/sum(bases[chrom%in%c("2L", "sim_2L")]),
                     coverage_across_reference_Dmel=mean(coverage_across_reference_Dmel),
                     totalreads=mean(totalreads),
                     PERCENT_DUPLICATION=mean(PERCENT_DUPLICATION)), list(sampleId)]



 melsim_hist <- ggplot(data=dt.sim, aes(mel_rate)) + geom_histogram()
 melsim <- ggplot(data=dt.sim, aes(y=coverage_across_reference_Dmel, x=totalreads*(1-PERCENT_DUPLICATION),
                  color=as.factor(mel_rate<.75))) +
 geom_point() +
 geom_text_repel(
   data=dt.sim[mel_rate<.75],
   aes(label = sampleId),
   size = 3,
   min.segment.length = 0,
   seed = 42,
   box.padding = 0.5,
   max.overlaps = Inf,
   arrow = arrow(length = unit(0.010, "npc")),
   nudge_x = .15,
   nudge_y = .5,
   color = "grey50"
 )

 melsim_hist + melsim
