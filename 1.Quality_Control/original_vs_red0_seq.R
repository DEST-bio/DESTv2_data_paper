#wc -l  /project/berglandlab/DEST/raw_reads/DrosEU_3_Jan2023/*.gz > reads.orig.num.txt
#wc -l  /project/berglandlab/DEST2.0_reads_deprecated/Feb2023_read_dump_ADDITIONAL/*.gz > reads.redo.num.txt

###
library(tidyverse)
library(vroom)
library(reshape2)
library(magrittr)

orig <- vroom("reads.orig.num.txt", col_names = F)
names(orig) = c("orig","ReadID")
redo <- vroom("reads.redo.num.txt", col_names = F)
names(redo) = c("redo","ReadID")

full_join(orig, redo) -> ori.redo

mean(ori.redo$orig/4, na.rm = T) -> glob_mean

ori.redo %>% 
  filter(!is.na(redo)) %>%
  melt(id = "ReadID") %>%
  ggplot(
    aes(
      x=ReadID,
      y=value/4,
      fill=variable
    )
  ) + coord_flip() +
  geom_bar(stat = "identity", 
           #position = "dodge"
           ) +
  geom_hline(yintercept = glob_mean)->
  pepi_redo

ggsave(pepi_redo, file = "pepi_redo.pdf")

