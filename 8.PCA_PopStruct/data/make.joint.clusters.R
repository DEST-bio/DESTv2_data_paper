### Make joint clusters
### 

library(tidyverse)
library(data.table)

meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"
samps <- fread(meta_git)
setDT(samps)


c4 <- get(load("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/8.PCA_PopStruct/data/nclust.4.sampleId.cluster.Rdata"))
names(c4)[2] = "cluster2.0_k4"
c5 <- get(load("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/8.PCA_PopStruct/data/nclust.5.sampleId.cluster.Rdata"))
names(c5)[2] = "cluster2.0_k5"
c8 <- get(load("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/8.PCA_PopStruct/data/nclust.8.sampleId.cluster.Rdata"))
names(c8)[2] = "cluster2.0_k8"

####
samps %>%
  full_join(c4) %>%
  full_join(c5) %>%
  full_join(c8) ->
  samps.clust

names(samps.clust)[which(names(samps.clust) == "cluster")] = "cluster1.0"

write.table(samps.clust, file = "dest_v2.samps_8Jun2023.csv", 
            append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
