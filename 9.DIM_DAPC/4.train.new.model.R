### Train new DAPC model
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)
library(data.table)
library(FactoMineR)

setwd("/scratch/yey2sn/DEST2_analysis/dapc_dims")

#### load PCA obj
pc.objc = get(load("../pca/pca.object.Rdata"))

pc.objc$eig[1:40,] %>%
  as.data.frame() %>%
  mutate(adj.perc = `percentage of variance`/40) %>%
  mutate(selected_snps = round(adj.perc*40000,0)) ->
  object.to.select.gims

object.to.select.gims$selected_snps %>% sum()
#### load and select SNPS

root = "/scratch/yey2sn/DEST2_analysis/pca/PCA_dimdescs"

GIMs.DEST2.0 =
foreach(i = 1:40, .combine = "rbind")%do%{
  message(i)
  
  fil <- paste("PC.",i,".new.mod.correlations.txt", sep = "")
  
  fil.root = paste(root, fil, sep = "/")
  
  tmp <- fread(fil.root) %>%
          arrange(V2) 
  
  tmp[1:object.to.select.gims$selected_snps[i],] ->
    tmp.sel
  
  names(tmp.sel) = c("cor","p.val.cor", "SNP_id", "PC")
  
  return(tmp.sel)
}

####
save(GIMs.DEST2.0, file = "DEST.2.0.GIMS.Rdata")

###
message("now loading AF")
AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))


####
grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps.o.nsim.noBeij
#grep( "SIM" , dat.o.nsim.noBeij$sampleId)

AF.d[-c(sim.pos, Beijing.pos), GIMs.DEST2.0$SNP_id] -> gims.AF

gims.AF %>% 
  PCA(graph = F, ncp = 5) ->
  pca.object.gims

save(pca.object.gims, file = "pca.object.gims.Rdata")





