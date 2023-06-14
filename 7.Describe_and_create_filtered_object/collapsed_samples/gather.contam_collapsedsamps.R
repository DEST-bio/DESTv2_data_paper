#### ---> Re run the contamnation assesment for DEST 2.0
#### Collapsed samples

### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(Rsamtools)
library(tidyverse)
library(vroom)

fils <- system("ls collapsed_sampsSimContam.*", intern = T)

collapse_contam = foreach(fi = fils, .combine = "rbind")%do%{
  
  tmp <- get(load(fi))
  
  return(tmp)
  
}

collapse_contam %>%
  group_by(sampleId) %>%
  summarise(MAPPED_eff = mean(mappingEffort, na.rm = T),
            SimCont.Norm = mean(propSimNorm, na.rm = T) ) ->
  contams.F

save(contams.F, file = "contams.F.collect.Rdata")

