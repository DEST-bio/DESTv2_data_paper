### libraries
  library(data.table)
  library(ggplot2)

### load samps
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_8Jun2023.csv")

### seasonal set
  load("~/seasonalpair.pca.meta.Rdata")
  ss_noCore20 <- seasonalpair.pca.meta[delta.T.mag=="Steep"][delta.T.sign==-1][Core20_sat==F]$sampleId
  ss_Core20 <- seasonalpair.pca.meta[Core20_sat==T]$sampleId

### make plot

  ggplot(data=samps[year>2009]) +
  geom_point(aes(y= lat,x= jday, group=city), size = 1.2) +
  geom_point(data=samps[year>2009 & sampleId%in%ss_noCore20], aes(y= lat,x= jday, group=city), size = 1.5, color="orange") +
  geom_point(data=samps[year>2009 & sampleId%in%ss_Core20], aes(y= lat,x= jday, group=city), size = 1.5, color="blue", symbol="diamond") +
  ylim(29,65) +
  xlim(100,350) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
  facet_grid(~year) + xlab("Julian Day") + ylab("Latitude")
