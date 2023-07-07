##### This is an R-script
##### Jcbn, June 8, 2023

### --> working forlder: /scratch/yey2sn/DEST2_analysis/admix_samps

##### Visualize F3 results

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

#### Obtain results

root_folder <- "/scratch/yey2sn/DEST2_analysis/admix_samps/Jun8.admix"
files_r <- system(paste("ls",root_folder, sep = " "), inter = T)

# Now load all results using a foreach loop with rbind as .combine
f3_results <- foreach(i=files_r, .combine = "rbind")%do%{

message(i)
tmp <- get(load(paste(root_folder,i, sep = "/")))
return(tmp)

}
setDT(f3_results)

# Now bring in metadata
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"
samps <- fread(meta_git)
setDT(samps)

samps %<>% mutate(focal = sampleId)
left_join(f3_results, samps) -> f3_meta

save(f3_meta, file = "f3_meta.joint.R")