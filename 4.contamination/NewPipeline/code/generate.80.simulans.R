##### Create SRA sampler for melanogaster

library(vroom)
library(tidyverse)

####
setwd("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/contamination/NewPipeline/")

sim.all.sea <- vroom("Simulans.SRAruns.txt")

sim.all.sea %>%
group_by(Strain) %>%
  slice_head() ->
  sim.subsam

#####
sim.subsam[sample(dim(sim.subsam)[1], 80),] -> 
  samps.80

samps.80 %>%
  dplyr::select(Run) ->
  o.file

#######
write.table(o.file, file = "Dsim.80.txt",
            col.names = T, row.names = F,
            quote = F)



