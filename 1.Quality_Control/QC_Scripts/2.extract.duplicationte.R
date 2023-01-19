### Mine the duplication rates
### 

library(vroom)
library(foreach)
library(tidyverse)
library(magrittr)






qcs <- vroom("qcs.batch.txt", col_names = F)

dup.rat.head <- read.table( pipe( paste("sed -n 7p ", "/project/berglandlab/DEST/dest_mapped/DESTv2/" , 
                                   paste(qcs$X1[1], "/", sep = ""),
                                   paste(qcs$X1[1], ".mark_duplicates_report.txt", sep = ""), 
                                   sep = "")))

dup.rate <- foreach(i=1:dim(qcs)[1], .combine = "rbind")%do%{
  
dup.rat <- read.table( pipe( paste("sed -n 8p ", "/project/berglandlab/DEST/dest_mapped/DESTv2/" , 
                    paste(qcs$X1[i], "/", sep = ""),
              paste(qcs$X1[i], ".mark_duplicates_report.txt", sep = ""), 
              sep = ""))
)

dup.rat %<>% as.data.frame() %>% mutate(samp = qcs$X1[i])
return(dup.rat)
}

names(dup.rate)[1:10] = dup.rat.head

names(qcs) = c("samp","path", "batch")

dup.rate %>%
  left_join(qcs) -> data.pcr.dups

data.pcr.dups %>%
  ggplot(aes(x=batch, y=PERCENT_DUPLICATION, fill = batch)) + 
  geom_boxplot() ->
  PERCENT_DUPLICATION
ggsave(PERCENT_DUPLICATION, file = "PERCENT_DUPLICATION.png", w = 5, h =4)

t.test(PERCENT_DUPLICATION~batch, data=data.pcr.dups)

