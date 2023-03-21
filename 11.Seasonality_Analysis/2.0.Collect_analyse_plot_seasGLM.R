#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
library(santoku)

#####

files.seas = system("ls ./GLM_out", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_out"

seas.p = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/"))) %>%
      mutate(perc_dat = nObs_i/nObs_tot)

    #perform binning with custom breaks
    tmp %>% 
      filter(perc_dat > 0.85) %>%
      mutate(AF_bin = chop(p_lrt, breaks=c(
        seq(from=0, to = 0.009, by = 0.001),
        seq(from=0.01, to = 0.09, by = 0.01),
        seq(from=0.1, to = 1.0, by = 0.1))
      )) -> tmp.i
    
    tmp.i %>%
      group_by(AF_bin, perm, chr) %>%
      summarize(Nsp = n()) %>%
      mutate(fil = fil)->
    tmp.ag
    
    return(tmp.ag)
    
  }


seas.p %>%
  group_by(AF_bin,perm,chr) %>%
  summarize(Nsps.ag = sum(Nsp)) ->
  seas.p.ag

save(seas.p.ag, file = "seas.p.ag.rdata")


