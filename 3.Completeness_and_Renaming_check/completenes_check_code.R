#### Renaming scheme:
#
#

library(vroom)
library(tidyverse)

####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_14Feb2023.csv")

####
master.samps = vroom("dest_v2.samps_14Feb2023.csv", delim = ",")

####
## ===> Get samples currently mapped: 
mapped.files =
system("ls /project/berglandlab/DEST/dest_mapped/*/*/*.sync.gz | grep -v 'SNAPE' | grep 'masked' ", inter = T)

mapped.files %>%
  data.frame(fi=.) %>%
  separate(fi, into = c("emp","prj","lab","proj","map","version","sampleId","sync"), sep = "/") ->
  mapped.files.dag

####
master.samps %>%
  mutate(validated = case_when(sampleId %in% mapped.files.dag$sampleId ~ "validated",
                               TRUE ~ "pending")) %>%
  filter(validated == "pending") %>%
  dplyr::select(sampleId, sampleId_orig, set, library_result) %>% 
  as.data.frame()
