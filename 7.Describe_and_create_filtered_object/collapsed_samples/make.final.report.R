
### Make final report
### 
library(tidyverse)
library(reshape2)
#clust <- get(load("sampleId.cluster.Rdata"))

contam <- get(load("contams.F.collect.Rdata"))

covs <- get(load("collapsed.covs.miss.pcr.Rdata"))
covs %>%
  dcast(sampleId~Var, value.var = "Value") -> cast.covs

left_join(cast.covs,contam) -> QC_data

save(QC_data, file = "QC_data.collapse.Rdata")
