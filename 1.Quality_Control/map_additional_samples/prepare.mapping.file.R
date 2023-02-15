#### Prepare the sampling scheme
#### 

library(vroom)
library(tidyverse)
library(magrittr)

# ====> download the metadata
#
#

system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_14Feb2023.csv")

####
master.samps = vroom("dest_v2.samps_14Feb2023.csv", delim = ",")

### Now load in the sampled that have been sequenced multiple times to increase coverage
### 

add.samps = vroom("/project/berglandlab/DEST2.0_reads/Feb2023_read_dump_ADDITIONAL/Feb2023_download_master.txt")
add.samps %>%
  separate(File, into = c("read_Dir","file","compression"), sep = "\\.") %>%
  separate(read_Dir, into = c("SequencingId", "PE_dir"), sep = "_") ->
  add.samps.pro

add.samps.pro$SequencingId %>% unique -> samps.to.map

####
master.samps %>%
  filter(SequencingId %in% samps.to.map) %>% 
  dplyr::select(sampleId, SequencingId, nFlies) ->
  metadat.to.map

write.table(metadat.to.map, file = "metadat.to.map.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


