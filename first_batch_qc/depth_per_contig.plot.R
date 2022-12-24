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
 geom_point()

 melsim_hist + melsim
