#### Quantify numbers in DEST

library(tidyverse)
library(data.table)
library(magrittr)

#####
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

####
samps$library_type %>% table
samps$country %>% table

####
samps %>%
  filter(library_type %in% c("pooled", "Pooled")) %>% 
  filter(is.na(collapsedSamples)) %>%
  group_by(locality, year) %>%
  summarise(N=n()) %>% .$N %>% table -> nomnum

samps %>%
  filter(library_type %in% c("pooled", "Pooled")) %>% 
  filter(is.na(collapsedSamples)) %>%
  group_by(locality, year) %>%
  summarise(N=n()) %>% .$N %>% table %>% prop.table -> freqnum

data.frame(nomnum, freqnum)

nomnum[-c(1:2)] %>% sum
freqnum[-c(1:2)] %>% sum

####
samps$Recommendation %>% table




