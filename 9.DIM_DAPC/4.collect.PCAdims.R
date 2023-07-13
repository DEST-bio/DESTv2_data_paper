### Train new DAPC model
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

##setwd("/scratch/yey2sn/DEST2_analysis/dapc_dims")

#### load PCA obj
pc.objc = get(load("/scratch/yey2sn/DEST2_analysis/pca/pca.object.noInvnoCd.250thinned.Rdata"))

pc.objc$eig[1:40,] %>%
  as.data.frame() %>%
  mutate(adj.perc = `percentage of variance`/40) %>%
  mutate(selected_snps = round(adj.perc*40000,0)) ->
  object.to.select.gims

object.to.select.gims$selected_snps %>% sum()
#### load and select SNPS

root = "/scratch/yey2sn/DEST2_analysis/dapc_dims/pca_mining"

GIMs.DEST2.0 =
foreach(i = 1:40, 
.combine = "rbind")%do%{

  message(i)
  
  fil <- paste("PC.",i,".new.mod.noInvnoCd.250thinned.correlations.txt", sep = "")
  
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



