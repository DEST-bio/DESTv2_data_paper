### Collect simulans data
library(foreach)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)
library(readxl)

###

files = system("ls | grep 'Rdata'", intern = T)

sim.cont = 
foreach(i=1:length(files), .combine = "rbind")%do%{
  
  message(i)
  tmp <- get(load(files[i]))
  tmp %>% 
    filter(chr %in% c("2L","2R","3L", "3R")) %>%
    group_by(sampleId) %>%
    summarize(m.sim = mean(propSim)) -> rt
  
  return(rt)
  
}

save(sim.cont, file = "sim.contam.joint.Rdata")