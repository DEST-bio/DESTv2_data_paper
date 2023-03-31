library(ggplot2)

meta <- read.csv("/media/DEST2_NHM/data/dest_v2.samps_25Feb2023.csv", row.names=1)
nhm_inversion <- read.delim("/media/DEST2_NHM/results/PoolSNP_nhm_inversion.af", row.names=1)
for (col in colnames(nhm_inversion)){
  nhm_inversion[[col]] <- as.numeric(paste(nhm_inversion[[col]]))
}  

joined_table <- merge(meta[3:5], nhm_inversion, by = "row.names")
rownames(joined_table) <- rownames(meta)
data <- merge(meta[13:14], joined_table, by ="row.names")
data2 <- na.omit(data) 

inversions=c("In.3R.Payne", "In.2L.t", "In.2R.Ns", "In.3R.C", "In.3R.K", "In.3R.Mo", "In.3L.P")

##global for all per continent (latitude)
#par(mfrow=c(2,3))
for (i in inversions){
  print(i)
  pdf(paste("/media/DEST2_NHM/results/Global_Latitude_",i,".pdf", sep=""),width=12,height=8)
  print(ggplot(data2,aes(x=lat,y=data2[[i]],col=continent))+
    geom_point()+
    geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ x, alpha=.15, aes(fill=continent)) +
    guides(fill=FALSE) + ggtitle(paste("Continents - ",i, sep=""))+theme(plot.title = element_text(hjust=0.5, size=24)) + theme(axis.text=element_text(size=24), axis.title=element_text(size=24))+  ylab(paste(i)))
  dev.off()
  }



##global for all per continent (longitude)
pdf("/media/DEST2_NHM/results/Global_Longitude.pdf",width=12,height=8)
#par(mfrow=c(2,3))
for (i in inversions){
  print(i)
  pdf(paste("/media/DEST2_NHM/results/Global_Longitude_",i,".pdf", sep=""),width=12,height=8)
  print(ggplot(data2,aes(x=long,y=data2[[i]],col=continent))+
          geom_point()+
          geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ x, alpha=.15, aes(fill=continent)) +
          guides(fill=FALSE) + ggtitle(paste("Continents - ",i, sep=""))+theme(plot.title = element_text(hjust=0.5)) +  ylab(paste(i)))
  dev.off
  }


## Additional Analysis for North American and European samples
NAm <- subset(data, continent=="North_America" )
Europe <- subset(data, continent=="Europe" )

# North America per Year
N12t18 <- subset(NAm, year>=2012 & year <=2018)
N12t18$year <- as.factor(paste(N12t18$year))
ggplot(N12t18,aes(x=lat,y=In.3R.Payne,col=year))+ geom_point()+ geom_smooth(method = "glm", method.args = list(family="binomial"), formula = y ~ x, alpha=.15, aes(fill=year))+  guides(fill=FALSE) + ggtitle("North American Samples: 2012-2018")+theme(plot.title = element_text(hjust=0.5)) 

E14t21 <- subset(Europe, year>=2014 & year <=2021)
E14t21$year <- as.factor(paste(E14t21$year))
ggplot(E14t21,aes(x=lat,y=In.3R.Payne,col=year))+ geom_point()+ geom_smooth(method = "glm", method.args = list(family="binomial"),alpha=.15, aes(fill=year))+  guides(fill=FALSE) + ggtitle("European Samples: 2014-2021")+theme(plot.title = element_text(hjust=0.5)) 

