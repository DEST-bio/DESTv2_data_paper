### libraries
  library(data.table)
  library(readxl)
  library(ggplot2)

### load ref genome based estimate
  load("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/1.Quality_Control/reseq/reseq.Rdata")
  setnames(m, "sim_alan", "Dsim_readDepth")

### load k-mer estimate
  kmer <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/contamination/dest_kmer_screening.xls"))
  kmerl <- melt(kmer, id.vars="Sample")

  mk <- merge(m, kmer[,c("Sample", "Dmela", "Dsimu"), with=F], by.x="sampleId", by.y="Sample", all=T)

### save
  save(mk, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/contamination/kmer_readDepth_data.Rdata")

  contam <- ggplot(data=mk, aes(x=(Dsim_readDepth), y=(Dsimu/100))) + geom_point(color="white", fill="black", shape=21, alpha=.9, size=3) +
  geom_abline(slope=1, intercept=0) + scale_x_log10() + scale_y_log10() +
  ylab("log10(D.sim contam): kmer") + xlab("log10(D.sim contam): read depth") + annotation_logticks() +
  geom_text_repel(data=mk[Dsim_readDepth>.02 & (Dsimu/100)<.01], aes(label=sampleId), size=2)

  ggsave(contam, file="~/contam.pdf")


  ggplot(data=mk[(Dsimu*nFlies/100)>1], aes(x=Dsimu*nFlies/100)) + geom_histogram(binwidth=.15)
