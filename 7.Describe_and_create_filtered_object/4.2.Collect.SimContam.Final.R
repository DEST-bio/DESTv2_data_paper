#### ---> Re run the contamnation assesment for DEST 2.0
#### 

### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(Rsamtools)
library(tidyverse)
library(vroom)

###
samps <- vroom("dest_v2.samps_25Feb2023.csv")
setDT(samps)
samps = samps[set!="dgn"]

####
fns <- system("ls /scratch/yey2sn/DEST2_analysis/filtering/sim_contam_final/*.Rdata", intern=T)

contams = 
foreach(i=fns, .combine = "rbind")%do%{
  tmp <- get(load(i))
}

contams %>%
  group_by(sampleId) %>%
  summarise(MAPPED_eff = mean(mappingEffort, na.rm = T),
            SimCont.Norm = mean(propSimNorm, na.rm = T) ) ->
  contams.F

save(contams.F, file = "Mean.contamination.Final.Rdata")

