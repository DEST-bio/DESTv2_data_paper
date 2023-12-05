library(ggplot2)
library(data.table)
samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

seasonal.sets <- samps[Recommendation=="Pass"][continent%in%c("Europe", "Asia")][set!="DrosRTEC"][,list(.N,
                                                               season=c("spring", "fall"),
                                                               sampleId=c(sampleId[which.min(jday)], sampleId[which.max(jday)]),
                                                               jday=c(jday[which.min(jday)], jday[which.max(jday)]),
                                                               diff=c(jday[which.min(jday)]-jday[which.max(jday)]),
                                                               lat=lat[which.min(jday)]),
                                                            list(country, locality, year, loc_rep, cluster2.0_k5)][N>1][abs(diff)>30]

collection_fig <- ggplot(data=seasonal.sets, aes(x=jday, y=lat, group=interaction(locality, year))) +
geom_line() + ylab("Latitude") + xlab("Julian Day")
