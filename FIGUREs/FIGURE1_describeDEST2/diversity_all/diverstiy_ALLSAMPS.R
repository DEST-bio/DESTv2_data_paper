### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)
  registerDoMC(4)
  library(ggbeeswarm)

### load npstat output
  dat <- readRDS("/Users/alanbergland/all.npstat.res.samplename.rds")

### load samps
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_8Jun2023.csv")
  setnames(samps, "cluster2.0_k4", "cluster")

### cycle and collect; use largest window size
  div <- foreach(i=1:length(dat[[1]][[3]]), .combine="rbind", .errorhandling="remove")%dopar%{
      #i <- 89
      tmp <- as.data.table(dat[[1]][[3]][i])
      setnames(tmp, "pool.lst.k..1.", "sampleId")
      tmp[,list(rd=mean(read_depth), W=mean(Watterson), Pi=median(Pi)), list(auto=as.factor(CHR!="X"), sampleId)]
  }

### merge

  div <- merge(div[auto==T], samps[,c("sampleId", "sampleId_orig", "continent", "country", "Recommendation", "cluster", "lat", "jday", "year"), with=F], by="sampleId", all.y=T)

  div[,cluster_afr:=as.character(cluster)]

  div[cluster_afr=="3", cluster_afr:="North Africa"]
  div[cluster_afr=="1", cluster_afr:="Sub-Saharan Africa"]
  div[!grepl("Africa", cluster_afr), cluster_afr:=NA]

### basic plot
  ggplot(data=div[Recommendation=="Pass"][auto==T], aes(x=interaction(continent), y=Pi)) + geom_jitter(width=.05)

### load in old data
  old_div <- fread("/Users/alanbergland/Documents/GitHub/data-paper/Figure6/Figure6_data.txt")[FILE=="PoolSNP"]
  divm <- merge(div[Recommendation=="Pass"], old_div, by.x="sampleId_orig", by.y="POP", all.x=T)
  divm[is.na(Pi.x)]

  divm[,Pi:=Pi.x]
  divm[is.na(Pi.x), Pi:=Pi.y]

## tidy
  divm <- divm[!is.na(continent)]
  divm[cluster=="1", continent:=cluster_afr]
  divm[cluster=="3" & continent=="Africa", continent:=cluster_afr]
  divm[country=="China", continent:="East Asia"]
  divm[continent=="Asia", continent:="West Asia"]
  divm[continent=="Europe" & cluster=="3", continent:="West Europe"]
  divm[continent=="Europe" & cluster=="2", continent:="East Europe"]

  divm[,continentf:=factor(continent, levels=c("Sub-Saharan Africa", "North Africa", "West Europe", "East Europe", "West Asia", "East Asia", "North_America", "South_America", "Oceania"))]
  divm[is.na(continentf)]

### plot
  save(divm, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE1_describeDEST2/diversity_all/diversity_dest2.Rdata")

  div_plot <- ggplot(data=divm, aes(x=(continentf), y=Pi, color=continentf)) +
  geom_jitter(width=.15, alpha=.75) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
  xlab("") + ylab("Theta Pi")
  ggsave(div_plot, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE1_describeDEST2/diversity_all/diversity_all.pdf")


  ggplot(data=divm[continent=="North_America"][grepl("US_Vir_Cha", sampleId)], aes(x=jday, y=Pi, color=continentf)) + geom_point() + facet_grid(~continentf)

  summary(lm(Pi~lat, divm[continent=="West Europe"]))
  summary(lm(Pi~lat, divm[continent=="East Europe"]))
  summary(lm(Pi~lat, divm[continent=="North_America"]))
  summary(lm(Pi~jday, divm[continent=="North_America"][grepl("US_Vir_Cha", sampleId)]))

  divm[is.na(Pi.x)]
