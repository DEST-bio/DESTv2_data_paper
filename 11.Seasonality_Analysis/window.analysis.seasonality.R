#### Collect and analyse the seasonality model -- in widnows

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)

###

glm.out <- get(load("seas.p.glm.2L.Rdata"))

# generate a master index for window analysis
### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(glm.out, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L"
                        #,"2R", "3L", "3R"
                        ),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- glm.out %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
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

    
  }

