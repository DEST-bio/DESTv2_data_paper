### Characterize Dataset
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

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

#####
#####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- vroom("dest_v2.samps_25Feb2023.csv")

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

####
#### ----> Generate Basic Stats
####
####
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))
filtering.dt %<>%
  filter(is.na(libs))


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
####
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
dat <- ad$data/dp
dim(dat)
colnames(dp) <- paste(seqGetData(genofile, "chromosome"),
                       seqGetData(genofile, "position") 
                       , sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

#### Estimate NA %
dim(dp)
rowSums(is.na(dp))/4570964 -> Missing.data.calc
rowMeans(dp, na.rm = T) -> Means.Cov

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


#