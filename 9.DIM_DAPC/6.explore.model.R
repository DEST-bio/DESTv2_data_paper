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


