##### Create SRA sampler for melanogaster

library(vroom)
library(tidyverse)

####
setwd("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/contamination/NewPipeline/")

mel.all.sea <- vroom("melanogaster.SRA.txt")

mel.all.sea %>%
group_by(Strain) %>%
  slice_head() ->
  dgrp.subsam

#####
dgrp.subsam[sample(dim(dgrp.subsam)[1], 80),] -> 
  samps.80

samps.80 %>%
  dplyr::select(Run) ->
  o.file

#######
write.table(o.file, file = "Dmel.80.txt",
            col.names = T, row.names = F,
            quote = F)



