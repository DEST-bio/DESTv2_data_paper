###
###

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(vroom)

###
dat.o <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DEST2.flt.meanAF0.01.miss0.01.Rdata"))

#####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- vroom("dest_v2.samps_25Feb2023.csv")

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps

dat.o[samps$sampleId,] -> dat.o.nsim.noBeij

#########
#########
#########

dat.o.nsim.noBeij %>%
  PCA(graph = F) ->
  pca.object

save(pca.object, file = "pca.object.Rdata")

