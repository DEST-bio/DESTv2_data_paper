###
###

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(vroom)
library(data.table)

library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(zoo)
library(adegenet)
library(reshape2)
library("geosphere")
library(gmodels)
library(ggforce)
library(forcats)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

####
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"

samps <- fread(meta_git)
setDT(samps)
samps.p = samps[set!="dgn"][Recommendation == "Pass"][collector != "Fournier-Level et al"]
samps.c = samps[!is.na(collapsedSamples)]
samps.af = samps[set=="dgn"][continent %in% c("Africa","North_America","Oceania","Europe") ]

samps = rbind(samps.p, samps.af, samps.c)

samps %>%
  group_by(set, continent) %>%
  summarise(N= n())

######
######
######
######

genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

####
set.seed(12345)
chr.vec = c("all","2L","2R","3L","3R")

pca.var.exp = list()

PCA.results.df =
foreach(i = 1:length(chr.vec),
        .combine = "rbind")%do%{
          
          seqResetFilter(genofile)
          
          chr.i = chr.vec[i]
          
          message(paste(chr.i, i, sep = "/"))
          
          if(chr.i == "all"){
            snps.dt.tmp = snps.dt
          } else if(chr.i != "all"){
            snps.dt.tmp = filter(snps.dt, chr == chr.i )
          }
               
  seqSetFilter(genofile, 
               sample.id=samps$sampleId, 
               variant.id=snps.dt.tmp[sample(dim(snps.dt)[1], 20000 )]$variant.id)
          
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  ad.matrix = ad$data
  dat <- ad.matrix/dp
  dim(dat)
  colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                         seqGetData(genofile, "position") 
                         , sep="_")
  rownames(dat) <- seqGetData(genofile, "sample.id")
  
  dat %>% 
    PCA(graph = F, ncp = 3) ->
    pca.object.plot
  
  pca.object.plot$ind$coord %>%
    as.data.frame() %>% 
    mutate(sampleId = rownames(.)) %>%
    left_join(samps) %>%
    filter(!is.na(continent)) %>%
    dplyr::select(Dim.1,Dim.2,Dim.3,sampleId,country,continent,lat,long)%>%
    mutate(case = chr.i)->
    pca.meta
  
  pca.var.exp[[i]] = pca.object.plot$eig %>% 
    as.data.frame() %>% mutate(PC = rownames(.)) %>%
    mutate(case = chr.i)
  
  return(pca.meta)
        } ### close for each
########
save(PCA.results.df, file = "PCA.results.df.Rdata")
save(pca.var.exp, file = "pca.var.exp.Rdata")


##### ---->>> Cluster Analysis
##### ---->>> Cluster Analysis
##### ---->>> Cluster Analysis
##### ---->>> Cluster Analysis
##### ---->>> Cluster Analysis

world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1
  ) + theme_classic() -> base_world


##### ---->>> Cluster Analysis
setDT(PCA.results.df)
pca.meta = PCA.results.df %>% filter(case == "all")
dat.cluster <- pca.meta[,c("Dim.1","Dim.2","Dim.3")]

pdf("clust.viz.ALL.pdf")
fviz_nbclust(dat.cluster, kmeans, nstart = 25,  method = "gap_stat", nboot = 100) 
dev.off()

foreach(i=c(4,5,8))%do%{
  
  set.seed(123)
  km.res <- kmeans(dat.cluster, i, nstart = 25)
  # 3. Visualize
  
  c.dat.data = cbind(pca.meta, cluster=km.res$cluster)[,c("sampleId","cluster")]
  c.dat.plo = cbind(pca.meta, cluster=km.res$cluster)
  
  save(c.dat.data, 
       file = paste("nclust",i,"sampleId.cluster.Rdata", sep = "."))
  
  if(i == 4){
    SET = "Set1"
  } else  if(i == 5){
    SET = "Dark2"
  } else if(i == 8){
    SET = "Set2"
  }
    
  base_world + 
    geom_point(
      data = c.dat.plo,
      aes(x=long,
          y=lat,
          fill = as.factor(cluster)), size = 2.3, shape = 21
    ) +
    scale_fill_brewer(palette = SET) -> cluster.plot
  
  ggsave(cluster.plot, file = paste("nclust",i,"cluster.plot.pdf", sep = ".")
         , w = 7, h = 3.5)
  
}

#### Zoom on Europe
#### Zoom on Europe
#### Zoom on Europe
#### Zoom on Europe
#### Zoom on Europe
world2 <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world2) +
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-12, 41.00), ylim = c(32.00, 63.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(data =c.dat.plo , 
             aes(x=as.numeric(long), 
                 y=as.numeric(lat), 
                 fill = as.factor(cluster)), size = 2.3, shape = 21)  +
  xlab("Lon") + ylab("Lat") + 
  ggtitle("Europe under K=8") + 
  theme(legend.position = "none") + 
  scale_shape_manual(values = c(21,22,23,24)) -> Suture_ZoneEU

ggsave(Suture_ZoneEU, file = "Suture_ZoneEU.pdf",  w = 7, h = 3.5)
######3
ggplot(data = world2) +
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-125.15, -55), ylim = c(-10, 50.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(data =c.dat.plo , 
             aes(x=as.numeric(long), 
                 y=as.numeric(lat), 
                 fill = as.factor(cluster)), size = 2.3, shape = 21)  +
  xlab("Lon") + ylab("Lat") + 
  ggtitle("North America under K=8") + 
  theme(legend.position = "none") + 
  scale_shape_manual(values = c(21,22,23,24)) -> Suture_ZoneAM

ggsave(Suture_ZoneAM, file = "Suture_ZoneAM.pdf",  w = 7, h = 3.5)


