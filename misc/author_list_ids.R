###authorship master
library(tidyverse)
library(reshape2)
library(magrittr)
library(data.table)
library(foreach)

#######
author <- fread("/Users/jcnunez/Library/CloudStorage/GoogleDrive-joaquin.c.b.nunez@gmail.com/My\ Drive/Nunez\ Lab/Papers_and_Readings/Nunez_Lab_Papers/In_Progress/DESTv2.0/authorship_master.upd.txt")

#### institutes
author$Institutions %>% unique() -> affiliations

affi_upd = c()
for(i in 1:length(affiliations)){
    message(i)
    tmp = affiliations[i]
    strsplit(tmp, "\\|")[[1]] -> split.vec
    affi_upd = append(affi_upd, split.vec)
}

data.frame(aff =unique(affi_upd), numaf = seq(from=1, to=length(unique(affi_upd)))) -> master_Aff

###

author_Afi_code = c()
for(i in 1:length(author$Name)){
  message(i)
  tmp = author[i,]
  strsplit(tmp$Institutions, "\\|")[[1]] -> split.inst
  affs = which(master_Aff$aff %in% split.inst)
  item = paste(c(tmp$Name, paste(affs, collapse = ",")), collapse  = " ")
  author_Afi_code = append(author_Afi_code, item) 
}

###

awks = c()
for(i in 1:length(author$Name)){
  message(i)
  tmp = author[i,]
  strsplit(tmp$Name, " ")[[1]] -> split.name
  item = paste(split.name[length(split.name)], "to" , tmp$ACKNOWLEDGMENTS)
  awks = append(awks, item) 
}

funds = c()
for(i in 1:length(author$Name)){
  message(i)
  tmp = author[i,]
  strsplit(tmp$Name, " ")[[1]] -> split.name
  item = paste(split.name[length(split.name)], "was supported by" , tmp$Funding)
  funds = append(funds, item) 
}