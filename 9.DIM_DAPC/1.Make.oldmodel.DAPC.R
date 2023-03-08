#### DAPC analysis
#### 

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

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

## 1. Get the data
system("wget https://github.com/DEST-bio/data-paper/blob/main/Figures7_and_9/data/AIM_SNPs.Rdata?raw=true -O AIM_SNPs.Rdata ")

## load data
load("AIM_SNPs.Rdata")

AIMS_Subset %>%
  separate(SNPid, into = c("chr", "pos", "snpid_dest1"), sep = "_") %>%
  mutate(DIM_mark = paste(chr, pos, sep = "_"))->
  AIMS_Subset.m

####
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))

#####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

samps <- vroom("dest_v2.samps_25Feb2023.csv")

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))
####
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
dat <- ad$data/dp
dim(dat)

#
colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                       seqGetData(genofile, "position") #, paste("snp", seqGetData(genofile,
                                                        #                          "variant.id"), sep = "")
                       , sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")
#

DEST1_markers = which(colnames(dat) %in% AIMS_Subset.m$DIM_mark)

dat[,DEST1_markers] -> dat_DEST1_markers
dat_DEST1_markers.naImp = na.aggregate(dat_DEST1_markers)

### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###
#dat_DEST1_markers.naImp[1:10,1:10]

samps %>%
  filter(sampleId %in% rownames(dat_DEST1_markers.naImp)) ->
  samps.m

samps.m$country = gsub("USA", "United States", samps.m$country)
samps.m$country = gsub("w501", "United States", samps.m$country)

samps.m %>%
  filter(set %in% c("DrosRTEC","DrosEU","dgn") ) ->
  samps.dest.samps.1

#### STATE MODEL BEGIN
#### STATE MODEL BEGIN
#### STATE MODEL BEGIN
#### STATE MODEL BEGIN
#### STATE MODEL BEGIN
#### STATE MODEL BEGIN

samps.dest.samps.1$province %>% 
  table %>%
  .[which(. > 1)] %>%
  names -> count_to_use

samps.dest.samps.1 %>%
  .[which(.$province %in% count_to_use),] -> samps_filt.d1

samples=dat_DEST1_markers.naImp %>%
  .[which(rownames(.) %in% samps_filt.d1$sampleId),]

grps=samps_filt.d1$province

####
DAPC_model.v1 <- xvalDapc(samples,
                       grps,
                       n.pca.max = 300,
                       training.set = 0.9,
                       result = "groupMean",
                       center = TRUE,
                       scale = FALSE,
                       n.pca = NULL,
                       n.rep = 30,
                       xval.plot = TRUE)

#save(DAPC_model.v1, file = "DAPC_model.v1.Rdata")
###
load("./DAPC_model.v1.Rdata")

samps.m %>%
  filter(set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ) ->
  samps.dest.samps.NewSamps

New_samples=dat_DEST1_markers.naImp %>%
  .[which(rownames(.) %in% samps.dest.samps.NewSamps$sampleId),]

predict.dapc(DAPC_model.v1$DAPC, newdata=New_samples) -> DEST1_predictions

###
samps$country = gsub("USA", "United States", samps$country)
samps$country = gsub("UK", "United Kingdom", samps$country)

###
samps %>%
  group_by(real.province=province) %>%
  summarize(r.long = mean(long),
            r.lat = mean(lat)
  ) -> real.lat.long

samps %>%
  group_by(pred.province=province) %>%
  summarize(p.long = mean(long),
            p.lat = mean(lat)
  ) -> pred.lat.long

samps %>%
  group_by(country, continent, real.province=province) %>%
  slice_head() -> countries


###

DEST1_predictions$posterior %>% 
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>% 
  left_join(samps.dest.samps.NewSamps[,c("sampleId","province")]) %>% 
  melt(id = c("sampleId","province")) ->
  all.preds.Dest1

save(all.preds.Dest1, file = "./all.preds.Dest1.Rdata")
load("all.preds.Dest1.Rdata")

### Top Predictions
all.preds.Dest1 %>%
  group_by(sampleId) %>%
  slice_max(value) %>% 
  dplyr::select(sampleId, 
                real.province = province, 
                pred.province = variable, 
                posterior = value ) %>% 
  left_join(real.lat.long) %>%
  left_join(pred.lat.long) %>% 
  group_by(sampleId) %>%
  mutate(hav_d = distHaversine(
    matrix(c(r.long, p.long, r.lat, p.lat),nrow = 2))/1000) %>% 
  left_join(dplyr::select(countries, country,continent,real.province), by = c("real.province")) ->
  hav.dist.obj.top.pred

save(hav.dist.obj.top.pred, file = "hav.dist.obj.top.pred.Rdata")
load("hav.dist.obj.top.pred.Rdata")

####### <<-=- STATE MODEL END
####### <<-=- STATE MODEL END
####### <<-=- STATE MODEL END
####### <<-=- STATE MODEL END
####### <<-=- STATE MODEL END
####### 
###
###### city model
samps.dest.samps.1$city %>% 
  table %>%
  .[which(. > 1)] %>%
  names -> city_to_use

samps.dest.samps.1 %>%
  .[which(.$city %in% city_to_use),] -> samps_filt.d1.cit

samples.c=dat_DEST1_markers.naImp %>%
  .[which(rownames(.) %in% samps_filt.d1.cit$sampleId),]

grps.c=samps_filt.d1.cit$city

####
DAPC_model.v1.city <- xvalDapc(samples.c,
                               grps.c,
                          n.pca.max = 300,
                          training.set = 0.9,
                          result = "groupMean",
                          center = TRUE,
                          scale = FALSE,
                          n.pca = NULL,
                          n.rep = 30,
                          xval.plot = TRUE)

samps.m %>%
  filter(set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ) ->
  samps.dest.samps.NewSamps

New_samples=dat_DEST1_markers.naImp %>%
  .[which(rownames(.) %in% samps.dest.samps.NewSamps$sampleId),]

predict.dapc(DAPC_model.v1.city$DAPC, newdata=New_samples) -> DEST1_predictions.city

###
samps$country = gsub("USA", "United States", samps$country)
samps$country = gsub("UK", "United Kingdom", samps$country)

###
samps %>%
  group_by(real.city=city) %>%
  summarize(r.long = mean(long),
            r.lat = mean(lat)
  ) -> real.lat.long.city

samps %>%
  group_by(pred.city=city) %>%
  summarize(p.long = mean(long),
            p.lat = mean(lat)
  ) -> pred.lat.long.city

samps %>%
  group_by(country, continent, real.city=city) %>%
  slice_head() -> countries.city


###

DEST1_predictions.city$posterior %>% 
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>% 
  left_join(samps.dest.samps.NewSamps[,c("sampleId","city")]) %>% 
  melt(id = c("sampleId","city")) ->
  all.preds.Dest1.city

save(all.preds.Dest1.city, file = "./all.preds.Dest1.city.Rdata")
load("all.preds.Dest1.city.Rdata")

### Top Predictions
all.preds.Dest1.city %>%
  group_by(sampleId) %>%
  slice_max(value) %>% 
  dplyr::select(sampleId, 
                real.city = city, 
                pred.city = variable, 
                posterior = value ) %>% 
  left_join(real.lat.long.city) %>%
  left_join(pred.lat.long.city) %>% 
  group_by(sampleId) %>% 
  mutate(hav_d = distHaversine(
    matrix(c(r.long, p.long, r.lat, p.lat),nrow = 2))/1000) %>% 
  left_join(dplyr::select(countries.city, country,continent,real.city), by = c("real.city")) ->
  hav.dist.obj.top.pred.city

save(hav.dist.obj.top.pred.city, file = "hav.dist.obj.top.pred.city.Rdata")
load("hav.dist.obj.top.pred.city.Rdata")

####### <<-=- CITY MODEL END
#######
###


#### 
#### 
#### #### PLOTs
#### PLOT
#### PLOT
#### PLOT


#Graph maps
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#### The case of Australia
#####
#### AUSTRALIA
filter(hav.dist.obj.top.pred.city, country == "Australia" ) %>%
  left_join(samps) -> aus.data

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-125.15, 149.99), ylim = c(-48.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_curve(aes(x = r.long, y = r.lat, xend = p.long, yend = p.lat, color = city),
             alpha = 0.7,
             data = aus.data ) +
  facet_wrap(~province) -> australia.plot
ggsave(australia.plot, file = "australia.plot.pdf", h = 4, w = 9)

#### SA
filter(hav.dist.obj.top.pred.city, continent == "South_America" ) %>%
  left_join(samps) -> SA.data

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-125.15, 149.99), ylim = c(-48.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_curve(aes(x = r.long, y = r.lat, xend = p.long, yend = p.lat, color = collector),
             alpha = 0.7,
             data = SA.data ) +
  facet_wrap(~country) -> SA.data.plot
ggsave(SA.data.plot, file = "SA.data.plot.pdf", h = 4, w = 9)


