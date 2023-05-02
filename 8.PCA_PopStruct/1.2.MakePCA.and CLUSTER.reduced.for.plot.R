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
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.qc_merge.csv"

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
seqSetFilter(genofile, sample.id=samps$sampleId, 
             variant.id=snps.dt[sample(dim(snps.dt)[1], 50000 )]$variant.id)

########
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
ad.matrix = ad$data
dat <- ad.matrix/dp
dim(dat)

colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                      seqGetData(genofile, "position") 
                      , sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

####
####
####

dat %>% 
  PCA(graph = F, ncp = 40) ->
  pca.object.plot

save(pca.object.plot, file = "pca.object.plot.Rdata")

#### #### #### #### #### #### #### #### #### #### #### #### 

pca.object.plot$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) %>%
  filter(!is.na(continent))->
  pca.meta

####
pca.meta %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = continent
  )) +
  geom_point(shape = 21, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() ->
  PCA12.lat

ggsave(PCA12.lat, file = "PCA12.lat.pdf", w = 5, h =4)

###
pca.meta %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.3,
    fill = continent
  )) +
  geom_point(shape = 21) +
  geom_point(shape = 21, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() ->
  PCA13.lat

ggsave(PCA13.lat, file = "PCA13.lat.pdf", w = 5, h =4)

###
###
###
###
###
###
###
###
cor.test(~lat + Dim.1, dat = pca.meta)
cor.test(~long + Dim.1, dat = pca.meta)

cor.test(~lat + Dim.2, dat = pca.meta)
cor.test(~long + Dim.2, dat = pca.meta)

cor.test(~lat + Dim.3, dat = pca.meta)
cor.test(~long + Dim.3, dat = pca.meta)

###
pca.meta %>%
  filter(continent %in% c("Europe","North_America","South_America","Oceania")) %>%
  dplyr::select(Dim.1,Dim.2,Dim.3, lat, long, continent) %>%
  reshape2::melt(id = c("lat", "long", "continent")) %>% 
  ggplot(aes(
    x=lat,
    y=value,
    color = variable
  )) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  xlab("Latitude") +
  ylab("PCA projections") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~continent, scale = "free_x") ->
  lat.detail
ggsave(lat.detail, file = "lat.detail.pdf", w = 8, h = 2.6)

pca.meta %>%
  filter(continent %in% c("Europe","North_America","South_America","Oceania")) %>%
  dplyr::select(Dim.1,Dim.2,Dim.3, lat, long, continent) %>%
  reshape2::melt(id = c("lat", "long", "continent")) %>% 
  ggplot(aes(
    x=long,
    y=value,
    color = variable
  )) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  xlab("Longitude") +
  ylab("PCA projections") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~continent, scale = "free_x") ->
  long.detail
ggsave(long.detail, file = "long.detail.pdf", w = 8, h = 2.6)

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
    color = "black", fill = "lightgray", size = 0.1
  ) + theme_classic() -> base_world


##### ---->>> Cluster Analysis

dat.cluster <- pca.meta[,c("Dim.1","Dim.2","Dim.3")]

set.seed(123)
km.res <- kmeans(dat.cluster, 5, nstart = 25)
# 3. Visualize

c.dat.data = cbind(pca.meta, cluster=km.res$cluster)[,c("sampleId","cluster")]
c.dat.plo = cbind(pca.meta, cluster=km.res$cluster)

save(c.dat.data, 
     file = "sampleId.cluster.Rdata")


base_world + 
  geom_point(
    data = c.dat.plo,
    aes(x=long,
        y=lat,
        fill = as.factor(cluster)), size = 1.5, shape = 21
  ) +
  scale_fill_brewer(palette = "Set1") -> cluster.plot

ggsave(cluster.plot, file = "cluster.plot.pdf", w = 7, h = 3.5)

