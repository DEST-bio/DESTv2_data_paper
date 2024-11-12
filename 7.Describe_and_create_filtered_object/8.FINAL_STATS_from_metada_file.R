library(data.table)
library(tidyverse)

data <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv"

dest_dat <- fread(data)

dest_dat %>% 
  filter(Recommendation %in% c("Pass") ) ->
  PASSsamps

PASSsamps$country %>% table %>% length()
PASSsamps$locality %>% table %>% length()

PASSsamps %>%
  group_by(locality, year) %>%
  slice_head() %>%
  summarise(N = n()) %>% 
  group_by(locality) %>%
  summarise(Nloc = sum(N)) %>%
  group_by(Nloc) %>%
  summarise(Ntot = n()) %>% 
  filter(Nloc > 2) %>% .$Ntot %>% sum



dest1$nFlies %>% sum

dest_dat %>% 
  filter(set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ) ->
  dest2
dest2$nFlies %>% sum

