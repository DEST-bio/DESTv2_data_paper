### Explore Model of DEST 2.0
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

####
####
####
#sftp://rivanna.hpc.virginia.edu/scratch/yey2sn/DEST2_analysis/dapc_dims/xval_province.DEST2.0.Rdata

model2.0 <- get(load("xval_province.DEST2.0.Rdata"))
model2.0$DAPC$pca.loadings %>% rownames() %>% sort -> GIMs
###
message("now loading AF")
AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))

####
#### Part 1. Split sample set..
#### 
i=1

which(colnames(AF.d) %in% GIMs) -> GIMS.id
AF.d[samps$sampleId[i],GIMS.id] %>%
  t() %>%
  as.data.frame() -> anchored.samp

anchored.samp %>%
  mutate(SNP_id = rownames(.)) %>%
  .[complete.cases(anchored.samp),] %>% 
  .$SNP_id ->
  GIMS.id.noNas

which(colnames(AF.d) %in% GIMS.id.noNas) -> GIMS.id.clean
AF.d[samps$sampleId[i],GIMS.id.clean] -> AF.d.anchor.clean
## Now get other samples
AF.d[samps$sampleId[-i],GIMS.id.clean] -> AF.d.OTHER.clean
AF.d.OTHER.clean_naImp = na.aggregate(AF.d.OTHER.clean)
####
rownames(AF.d.OTHER.clean_naImp) -> ids.used.in.modeltrain
samps %>%
  filter(sampleId %in% ids.used.in.modeltrain) ->
  samps.used.in.modeltrain

##########
##########
##########
dapc(AF.d.OTHER.clean_naImp, 
     grp=samps.used.in.modeltrain$province, 
     n.pca=100, n.da=75) ->
  model.lo_out





