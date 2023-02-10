### load in the datasets 
### 

library(vroom)
library(tidyverse)
library(foreach)

sim.samps <- vroom("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/contamination/NewPipeline/Guide_files/Dsim.80.txt")
mel.samps <- vroom("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/contamination/NewPipeline/Guide_files/Dmel.80.txt")

####
####

foreach(i = 1:50, .combine = "rbind")%do%{
  
  sim.share = 80*(i/50)
  mel.share = 80
   
}


