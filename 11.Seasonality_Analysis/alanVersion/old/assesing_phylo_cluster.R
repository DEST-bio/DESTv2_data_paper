### IDENTIFY PHYLO CLUSTER + PCA data
### 
### 

### libraries
library(data.table)
library(lubridate)
library(foreach)
library(SeqArray)
library(doMC)
registerDoMC(5)
library(tidyverse)
library(lme4)
library(FactoMineR)
library(factoextra)

### load data
### seasonal pairs
seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
### core20  
core.20 <- fread("./core20_samps.csv")
names(core.20)[1] = "sampleId_orig"
### dest samps  
samps <- fread("./dest_v2.samps_25Feb2023.csv")

core20.upd = left_join(core.20, samps[,c("sampleId", "sampleId_orig")])

###
### gds object
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds")

### sample metadata
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

### get basic index
data <- seqGetData(genofile, "annotation/info/AF")
seqResetFilter(genofile)
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAllele=seqGetData(genofile, "$num_allele"),
                     variant.id=seqGetData(genofile, "variant.id"))

snp.dt <- snp.dt[nAllele==2]
seqSetFilter(genofile, snp.dt$variant.id)
snp.dt[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]
##
snp.dt[global_af > 0.05][chr %in% c("2L","2R","3L","3R")] -> snp.dt.flt

snp.dt.subset[sample(dim(snp.dt.subset)[1],10000, replace = FALSE)] %>%
.[order(variant.id)] -> snp.dt.flt.samp

seqSetFilter(genofile, snp.dt.flt.samp$variant.id)
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

##### add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")
#### 
#### 
dat[which(rownames(dat) %in% seasonal.sets$sampleId),] -> dat.seas

####
####

dat.seas %>% 
  PCA(scale.unit = F, graph = F, ncp = 50) -> 
  seas.pop.PCAobj

seas.pop.PCAobj$ind$coord %>% 
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>% 
  left_join(., seasonal.sets) -> 
  PCA_coords_metadata

###
set.seed(123)
km.res <- kmeans(PCA_coords_metadata[1:50], 3, nstart = 25)
# 3. Visualize
data.frame(sampleId = rownames(dat.seas), cluster = km.res$cluster) -> phylocluster_data
save(phylocluster_data, file = "phylocluster_data.Rdata")

PCA_coords_metadata %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color=season,
    shape = as.factor(km.res$cluster)
  )) +
  geom_point(size = 1.9) ->
  seas.pca.pop

ggsave(seas.pca.pop, file = "seas.pca.pop.pdf")


#####
dat[which(rownames(dat) %in% core20.upd$sampleId),] -> dat.core20

dat.core20 %>% 
  PCA(scale.unit = F, graph = F, ncp = 50) -> 
  seas.core20.pop.PCAobj

seas.core20.pop.PCAobj$ind$coord %>% 
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>% 
  left_join(., seasonal.sets) -> 
  PCA_coords_metadata_Core20

PCA_coords_metadata_Core20 %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color=season,
    #shape = as.factor(km.res$cluster)
  )) +
  ggtitle("Core20") +
  geom_point(size = 1.9) ->
  seas.pca.pop.core20

ggsave(seas.pca.pop.core20, file = "seas.pca.pop.core20.pdf")
