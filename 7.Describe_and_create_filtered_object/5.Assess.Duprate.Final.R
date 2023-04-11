##### ---> Duplication rates
##### 

## libraries
library(vroom)
library(foreach)
library(tidyverse)
library(magrittr)
library(data.table)


samps <- vroom("dest_v2.samps_25Feb2023.csv")
setDT(samps)
samps = samps[set!="dgn"]

qcs <- system("ls /project/berglandlab/DEST/dest_mapped/*/*/*duplicates_report.txt", intern=T)

dup.rat.head <- read.table( pipe( paste("sed -n 7p ", 
                                        paste(qcs[1]) 
                                        )
                                  ))
####
dup.rate <- foreach(i=1:length(qcs), .combine = "rbind")%do%{
  
  
  
  dup.rat <- read.table( pipe( paste("sed -n 8p ", 
                                     paste(qcs[i]))
  ))
  
  dup.rat %<>% as.data.frame() %>% mutate(samp = str_split( qcs[i], "/")[[1]][7] )
  message(paste(i, dim(dup.rat)[2], sep = "-") )
  return(dup.rat)
}
names(dup.rate)[1:10] = dup.rat.head

save(dup.rate, file = "DuplicateRates.all.Rdata")

