#### Train new DIMS models
#### 

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(FactoMineR)
library(factoextra)
library(tidyverse)

print(args[1])

pc.objc = get(load("./pca.object.Rdata"))

temp_dimdesc = dimdesc(pc.objc, axes = as.numeric(args[1]), proba = 0.05)

temp_dimdesc[[c( paste("Dim",args[1], sep =".") )]] %>% .$quanti %>% as.data.frame() %>% mutate(snpID = rownames(.), PC = as.numeric(args[1])) -> corr_temp

write.table( corr_temp ,
             file = paste("PC",args[1],"new.mod.correlations.txt", sep = "."),
             sep = "\t",quote = F ,row.names = F, col.names = F, append = T)


