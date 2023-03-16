###
###

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(vroom)

###
dat.o <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
#####

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps.o.nsim.noBeij
#grep( "SIM" , dat.o.nsim.noBeij$sampleId)

### --- filter data ---

dat.o[samps.o.nsim.noBeij$sampleId, 
      #sample(dim(dat.o)[2], 10000) 
      ] -> dat.o.nsim.noBeij

###### Filter SNPs by thinning phisical linakge
source("/home/yey2sn/software/ThinLDinR_SNPtable.R")
# pickSNPs<-function(map, dist = 100000) ...
### --> Working on this <<<<< ----
colnames(dat.o) -> snps.ids

data.frame(SNP_id = snps.ids) %>%
  separate(SNP_id, into = c("chr","pos","id")) ->
  SNPs.df

SNPs.df$pos = as.numeric(SNPs.df$pos)

SNPs.df %>%
  arrange(chr, pos) ->
  SNPs.df.sort

SNPs.df.sort %>% head

### apply thinning
pickSNPs(SNPs.df.sort, dist = 5000 ) -> thinned.set
thinned.set %>% length()

SNPs.df.sort[thinned.set,] %>% 
  mutate(SNP_id = paste(chr, pos, id, sep = "_")) ->
  SNPs.df.sort.thinned
  
#########
#########
#########

which(colnames(dat.o.nsim.noBeij) %in%  SNPs.df.sort.thinned$SNP_id) -> SNP.selector

dat.o.nsim.noBeij[,SNP.selector] %>% 
  PCA(graph = F, ncp = 5) ->
  pca.object.plot

###save(pca.object, file = "pca.object.Rdata")

####

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
