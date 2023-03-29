#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
#####


#### Central files
files.seas = system("ls ./GLM_out.03.27.2023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_out.03.27.2023"

#### binning analysis

seas.p.bin = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/"))) %>%
      mutate(perc_dat = nObs_i/nObs_tot)

    #perform binning with custom breaks
    #tmp %>% 
    #  filter(perc_dat > 0.85) %>%
    #  mutate(AF_bin = chop(p_lrt, breaks=c(
    #  seq(from=0.001, to = 0.009, by = 0.001),
    #  seq(from=0.01, to = 0.09, by = 0.01),
    #  seq(from=0.1, to = 1.0, by = 0.1))
    #  )) -> tmp.i
    #
    #tmp.i %>%
    #  group_by(AF_bin, perm, chr) %>%
    #  summarize(Nsp = n()) %>%
    #  mutate(fil = fil)->
    #tmp.ag
    
    return(tmp)
    
  }

seas.p.bin %>%
  filter(perm == 0) %>%
  ggplot(aes(
    p_lrt
  )) +
  geom_histogram() ->
  hist.test

ggsave(hist.test, file = "hist.test.pdf")


seas.p.bin %>%
  group_by(AF_bin,perm,chr) %>%
  summarize(Nsps.ag = sum(Nsp)) ->
  seas.p.bin.ag

save(seas.p.bin.ag, file = "seas.p.bin.ag.Rdata")

#### --> collect object GLMS
#### 
seas.p.glm.2L = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/"))) %>%
      mutate(perc_dat = nObs_i/nObs_tot) %>%
      filter(chr == "2L")
    
    return(tmp)
  }

seas.p.glm.2L$pos = as.numeric(seas.p.glm.2L$pos)
seas.p.glm.2L %<>% arrange(chr, pos) 

save(seas.p.glm.2L, file = "seas.p.glm.2L.Rdata")


