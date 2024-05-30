#### f3 table

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)

#####
f3tab <- fread("/gpfs2/scratch/jcnunez/DEST2.0_analysis/f3_revist/f3_results.txt")

f3tab %<>%
mutate(signif = 
case_when( `Z-score` <  -1.65 ~ "S",
TRUE ~ "NS"))

f3tab %>% group_by(focal) %>%
count(signif) %>%
dcast(focal~signif) %>%
mutate(Psig = 100*(S/(NS+S) )) %>%
.[complete.cases(.),] -> proStest

f3tab %>% group_by(focal) %>%
slice_min(`Z-score`) %>%
dplyr::select(focal, european_parent, african_parent) -> parent

left_join(proStest, parent ) -> tablefinal