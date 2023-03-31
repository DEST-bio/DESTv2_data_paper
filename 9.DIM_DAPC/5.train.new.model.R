### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)
library(data.table)
library(FactoMineR)
library(adegenet)
library(zoo)

setwd("/scratch/yey2sn/DEST2_analysis/dapc_dims")

load("DEST.2.0.GIMS.Rdata")
#####
message("now loading AF")
AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
##
######
grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps.o.nsim.noBeij
#grep( "SIM" , dat.o.nsim.noBeij$sampleId)

AF.d[-c(sim.pos, Beijing.pos), GIMs.DEST2.0$SNP_id] -> gims.AF
##
##gims.AF %>% 
##  PCA(graph = F, ncp = 5) ->
##  pca.object.gims
##
##save(pca.object.gims, file = "pca.object.gims.Rdata")
### Remove singletons
samps.o.nsim.noBeij %>% 
  group_by(province) %>%
  summarize(N = n()) %>%
  filter(N == 1) -> singleton.cities

samps.o.nsim.noBeij %>%
  filter(!province %in% singleton.cities$province) ->
  gim.samps

gim.samps$province = gsub("Kiev City", "Kiev", gim.samps$province)

AF.d[gim.samps$sampleId, GIMs.DEST2.0$SNP_id] -> gims.flt.AF

### clean AF tables from NA
#apply(gims.flt.AF, 2, function(x) sum(is.na(x)) ) -> Count.of.NAs
message("now na.aggregate")

gims.flt.AF_naImp = na.aggregate(gims.flt.AF)


###
message("now training xval")
xval_province <- xvalDapc(gims.flt.AF_naImp, 
                          as.factor(gim.samps$province), 
                          n.pca.max = 300, 
                          training.set = 0.9,
                          result = "groupMean", 
                          center = TRUE, 
                          scale = FALSE, 
                          n.pca = NULL, 
                          n.rep = 30, 
                          xval.plot = FALSE)

save(xval_province, file = "xval_province.DEST2.0.Rdata")


