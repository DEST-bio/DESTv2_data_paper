#### Mine new DIMS models
#### 

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(FactoMineR)
library(factoextra)
library(tidyverse)

print(args[1])

pc.objc = get(load("/scratch/yey2sn/DEST2_analysis/pca/pca.object.noInvnoCd.250thinned.Rdata"))

#temp_dimdesc = dimdesc(pc.objc, axes = 1:40, proba = 0.05)

temp_dimdesc = dimdesc(pc.objc, axes = as.numeric(args[1]), proba = 0.05)

temp_dimdesc[[c( paste("Dim",args[1], sep =".") )]] %>% .$quanti %>% as.data.frame() %>% mutate(snpID = rownames(.), PC = as.numeric(args[1])) -> corr_temp

write.table( corr_temp ,
             file = paste("./pca_mining/PC",args[1],"new.mod.noInvnoCd.250thinned.correlations.txt", sep = "."),
             sep = "\t",quote = F ,row.names = F, col.names = F, append = T)

#save(corr_temp, file = "temp_dimdesc.noInvnoCd.250thinned.Rdata")

