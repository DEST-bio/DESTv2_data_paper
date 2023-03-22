### Characterize Data set
### 

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


#####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- fread("dest_v2.samps_25Feb2023.csv")


##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA
##### DESCRIBE ORIGINAL METADATA


samps %>% dim
samps %>% group_by(set) %>%
  summarize(N = n())

samps %>% group_by(continent) %>%
  summarize(N = n())
samps %>% group_by(country) %>%
  summarize(N = n())

samps %>% group_by(year) %>%
  summarize(N = n()) %>% tail


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#### ---> Geographical plot
ggplot(data = world) +
  geom_sf(fill= "grey100") +
  coord_sf(xlim = c(-127.15, 149.99), ylim = c(-48.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "lightblue")) +
  geom_point(data = samps,
             aes(
               x=long,
               y=lat,
               fill=set
             ), shape = 21, size = 2.8, color = "white") -> DEST2.0.map
ggsave(DEST2.0.map, file = "DEST2.0.map.pdf", h = 4, w = 9)

####
ggplot(data = world) +
  geom_sf(fill= "grey100") +
  coord_sf(xlim = c(-12.50, 41.00), ylim = c(33.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "lightblue")) +
  geom_point(data = samps,
             aes(
               x=long,
               y=lat,
               fill=set
             ), shape = 21, size = 2.8, color = "white") -> DEST2.0.map.EU
ggsave(DEST2.0.map.EU, file = "DEST2.0.map.EU.pdf", h = 4, w = 9)

##### --> Seasonal plot
samps %>% group_by(year, city)%>%
  summarise(N=n()) %>% .$N %>% table %>% prop.table()*100

samps %>% group_by(year, city)%>%
  summarise(N=n()) %>% arrange(-N) %>%
  filter(N > 1) -> multisamps

samps %>%
  filter(!(set %in% c("dgn"))) %>%
  filter(continent %in% c("Europe", "North_America", "Asia")) %>%
  mutate(Region = case_when(long < -20 ~ "1.Long < -20",
                            long > -20 & long < 20 ~ "2.20 < Long > -20",
                            long >  20 ~ "3.Long > -20",
  )) %>%
  filter(city %in% multisamps$city) %>%
  separate(sampleId, 
           into = c("CountryT", "RegionT", "CityT", "RepT", "CollecDate"),
           sep = "_") %>%
  ggplot(aes(
    x=as.Date(CollecDate),
    y=province,
    #color=country,
    group = paste(province)
  )) + geom_line() + geom_point() +
  facet_wrap(~Region, scales = "free" )->
  date.test.plot
ggsave(date.test.plot, file = "date.test.plot.pdf", w = 9, h = 5)

#####
#####
#####
#####
##### DESCRIBE ORIGINAL GDS
##### DESCRIBE ORIGINAL GDS
##### DESCRIBE ORIGINAL GDS
##### DESCRIBE ORIGINAL GDS
##### DESCRIBE ORIGINAL GDS

#### ----> Generate Basic Stats
####
####
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))
filtering.dt %<>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) 


#grep( "SIM" , samps$sampleId) -> sim.pos
#grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos
#samps[-c(sim.pos, Beijing.pos),] -> samps

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=filter(samps, set %in% c("dest_plus","DrosEU_3","DrosEU_3_sa"))$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))
snps.dt <- snps.dt[nAlleles==2][missing<.10]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %<>% filter(SNP_id %in% filtering.dt$SNP_id)
snps.dt %>% dim

####
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

########
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
ad.matrix = ad$data
dat <- ad.matrix/dp
dim(dat)

colnames(dp) <- paste(seqGetData(genofile, "chromosome"),
                      seqGetData(genofile, "position") 
                      , sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

#### Estimate NA %

rowSums(is.na(dp))/dim(dp)[2] -> Missing.data.calc
rowMeans(dp, na.rm = T) -> Means.Cov

####### COVERGAE DATA
Means.Cov.df =
  data.frame(Value = Means.Cov,
             Var = "Cov") %>% 
  mutate(sampleId = names(Means.Cov)) %>%
  left_join(samps)

### Missing data
Missing.data.df =
  data.frame(Value = Missing.data.calc,
             Var = "Miss") %>% 
  mutate(sampleId = names(Missing.data.calc)) %>%
  left_join(samps)
####

#save(Miss.Cov.joint, file = "Miss.and.Cov.joint.Rdata")

### PCR DUPLICATES
### 
 fns <- system("find /project/berglandlab/DEST/dest_mapped/ -name '*duplicates_report.txt'", intern=T)
 pcr <- foreach(fn=fns, .combine = "rbind")%dopar%{
   #fn <- fns[1]
   message(fn)
   tmp <- fread(fn, skip="LIBRARY", nrows=1)
   data.table(pcrDup=tmp$PERCENT_DUPLICATION, READ_PAIRS_EXAMINED=tmp$READ_PAIRS_EXAMINED, sampleId=tstrsplit(fn, "/")[[7]])
 }

  pcr %>%
   as.data.frame() %>%
   dplyr::select(Value = pcrDup, sampleId) %>% 
   mutate(Value=Value, Var = "PCRdup")  %>%
    left_join(samps) ->
    pcr.mod

#### CONTAMINATION   
sim.contam <- get(load("simulans_rates.Rdata")) ### found un github!!
sim.contam %>%
  dplyr::select(sampleId=samp, Value=simNorm ) %>% 
  mutate(Value=Value, Var = "SimCont") %>%
  left_join(samps) ->
  sim.contam.samps
  
####  
  rbind(Means.Cov.df, Missing.data.df, pcr.mod, sim.contam.samps) -> Miss.Cov.PCRdup.sim.joint
  
  save(Miss.Cov.PCRdup.sim.joint,
       file = "Miss.Cov.PCRdup.sim.joint.Rdata")
  
  ### Filters
  Miss.Cov.PCRdup.sim.joint %>%
    dplyr::select(Value,Var,sampleId) %>%
    dcast(sampleId~Var, value.var = "Value") %>% 
    mutate(keep = case_when(
      Cov > 10 & Miss < 0.1 & PCRdup < 1.0 & SimCont < 0.1 ~ "PASS",
      TRUE ~ "FAIL"
    )) %>% left_join(samps) ->
    keep.fail.samps
  
  save(keep.fail.samps,
       file = "keep.fail.samps.Rdata")
  ##### ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ### object to save ####
  
####  
  keep.fail.samps %>%
    ggplot(aes(
      x=keep,
      fill=set 
    )) +
    geom_bar(position = "dodge") ->
    keep.lost.plot
  ggsave(keep.lost.plot, file = "keep.lost.plot.pdf", w = 9, h = 6)
  
  
####
  Miss.Cov.PCRdup.sim.joint %>%
    filter(!is.na(set)) %>%
    ggplot(aes(
      x=Value,
      color=set
    )) +
    geom_histogram() +
    facet_wrap(~Var, scale = "free") ->
    hist.QC.stats
  ggsave(hist.QC.stats, file = "hist.QC.stats.pdf", w = 9, h = 6)
  
  #filter(Miss.Cov.PCRdup.sim.joint, Var == "Miss") %>% arrange(Value) %>% tail
  
####    
  
 ###
 data.frame(
   Met = c("median", "median","median" , "median", "tresh","tresh","tresh","tresh"),
   Var = c("Cov","Miss","PCRdup","SimCont","Cov","Miss","PCRdup","SimCont"),
   m = c(median(filter(Miss.Cov.PCRdup.sim.joint, Var == "Cov")$Value),
         median(filter(Miss.Cov.PCRdup.sim.joint, Var == "Miss")$Value),
         median(filter(Miss.Cov.PCRdup.sim.joint, Var == "PCRdup")$Value),
         median(filter(Miss.Cov.PCRdup.sim.joint, Var == "SimCont")$Value, na.rm = T),
         10, 0.1, 1.0, 0.1)
 ) -> m.dfs
 
 ####
 Miss.Cov.PCRdup.sim.joint %>%
   filter(Var == "Cov") %>%
   arrange(Value) %>%
   mutate(order.id = 1:737) %>%
   dplyr::select(sampleId, order.id) -> order.vec
##### 
 
 Miss.Cov.PCRdup.sim.joint %>%
   left_join(order.vec) %>%
   filter(!is.na(set)) %>%
   ggplot(aes(
     x=fct_reorder(sampleId, order.id),
     y=Value,
     color = set
   )) +
   geom_point() +
   geom_hline(data =m.dfs, aes(yintercept = m, linetype = Met ) ) +
   #geom_hline(yintercept = quantile(Means.Cov.df$Value)[2], linetype = "dashed" ) +
   #geom_hline(yintercept = quantile(Means.Cov.df$Value)[4], linetype = "dashed" ) +
   facet_grid(Var~set, scales = "free", space = "free") +
   scale_y_continuous(trans='log10') +
   theme_bw() +
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.pos = "bottom") ->
   means.plot
 ggsave(means.plot, file = "means.plot.pdf", w = 10, h = 4.0)
 
 
### ---> Ended here for latest figure! More granular figures below + the aggregation set
### 
### 
### 
### 
### 
### 
### 
### 
cbind(Means.Cov, Missing.data.calc) %>% 
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>% 
  left_join(samps) %>% 
  group_by(sampleId) %>%
  mutate(Neff = (Means.Cov * nFlies-1)/(Means.Cov +nFlies ) ) ->
  seq.stats.df
###
seq.stats.df %>%
  dplyr::select(sampleId, continent, country, collector, 
                Means.Cov, Missing.data.calc, set, nFlies, Neff) -> DEST2.0.stats.summary

save(DEST2.0.stats.summary,
     file = "DEST2.0.stats.summary.Rdata")

load("DEST2.0.stats.summary.Rdata")
#####
##### Some correlations
lm(data=filter(DEST2.0.stats.summary, 
               Missing.data.calc < 0.07),
   Missing.data.calc ~ Means.Cov
) %>% summary

cor.test( ~Missing.data.calc+ Means.Cov,
          data =filter(DEST2.0.stats.summary, Missing.data.calc < 0.07))


#### Some aggregations
DEST2.0.stats.summary %>%
  dplyr::select(!collector) %>%
  filter(set != "dgn") %>%
  melt(id = c("sampleId", "continent", "country", "set")) %>% 
  group_by(set,  variable) %>% 
  summarise(mean = median(value, na.rm = T)) %>%
  dcast(set ~ variable)

#####
DEST2.0.stats.summary %>%
  filter(set == "dest_plus") %>% 
  melt(id = c("sampleId", "continent", "country", "collector", "set")) %>% 
  group_by(collector,  variable) %>% 
  summarise(mean = ci((value), na.rm = T)[1],
            uci = ci((value), na.rm = T)[3],
            lci =  ci((value), na.rm = T)[2]
  ) -> DEST_plus_stats

DEST_plus_stats %>%
  ggplot(aes(
    x=collector ,
    y=mean,
    ymin = lci,
    ymax = uci
  )) +
  geom_errorbar(width = 0.5) +
  geom_point() +
  coord_flip() + 
  facet_grid(~variable, scales = "free") ->
  dest.plus.plot
ggsave(dest.plus.plot, file = "dest.plus.plot.pdf", h = 2.3, w = 8)


### plot
### plot
### plot
### plot
### plot
DEST2.0.stats.summary %>%
  filter(set != "dgn") %>%
  melt(id = c("sampleId", "continent", "country", "set")) %>% 
  ggplot(aes(
    x=set,
    y=value,
  )) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  coord_flip() +
  facet_wrap(~variable, scales = "free") -> stats.general
ggsave(stats.general, file ="stats.general.pdf", w = 6, h = 4)


##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET
##### DESCRIBE MERGE DATASET

#####
#####
#####
#####
#####

metadata <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))

### Load objects 
DPmatrix <- get(load("/project/berglandlab/DEST2.0_working_data/joint.dp.matrix.Rdata"))
ADmatrix <- get(load("/project/berglandlab/DEST2.0_working_data/joint.ad.matrix.Rdata"))
AFmatrix <- get(load("/project/berglandlab/DEST2.0_working_data/joint.AFs.Rdata"))

##### DESCRIBE DATASETS
#### Estimate NA %
dim(DPmatrix)[2]
rowSums(is.na(DPmatrix))/(dim(DPmatrix)[2]) -> Missing.data.calc
rowMeans(DPmatrix, na.rm = T) -> Means.Cov

####
cbind(Means.Cov, Missing.data.calc) %>% 
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>% 
  left_join(metadata) %>% 
  group_by(sampleId) %>%
  mutate(Neff = (Means.Cov * nFlies-1)/(Means.Cov +nFlies )) ->
  metadata.updated

metadata.updated %>%
  filter(collector == "Fournier-Level et al") %>%
  group_by(collector) %>%
   summarise(Neff = mean(Neff),
            mCov = mean(Means.Cov))

###
metadata.updated %>%
  filter(Missing.data.calc < 0.3) ->
  metadata.flt
metadata.flt %>% dim
#
metadata.updated %>%
  filter(Missing.data.calc > 0.3) %>% dim() 

metadata.flt %>%
  filter(sampleId %in% samps.to.retain) %>%
  group_by(set) %>%
  summarise(m = mean(Means.Cov),
            ne = mean(Neff),
            miss = mean(Missing.data.calc))

#### Apply filtering
samps.to.retain = metadata.flt$sampleId

###
DPmatrix[samps.to.retain,] -> DPmatrix.flt
ADmatrix[samps.to.retain,] -> ADmatrix.flt
AFmatrix[samps.to.retain,] -> AFmatrix.flt

###
save(DPmatrix.flt,
     file = "/project/berglandlab/DEST2.0_working_data/DPmatrix.flt.Rdata")
save(ADmatrix.flt,
     file = "/project/berglandlab/DEST2.0_working_data/ADmatrix.flt.Rdata")
save(AFmatrix.flt,
     file = "/project/berglandlab/DEST2.0_working_data/AFmatrix.flt.Rdata")
###

