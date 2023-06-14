#### ---> Create Filtering set information
#### 
library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr)
library(data.table)

### from github
qc.dat <- get(load("Miss.Cov.PCRdup.sim.joint.Rdata"))

qc.dat %>% 
  dcast(sampleId~Var, value.var = "Value") ->
  qc.dat.miss.cov.etc


####


## Set global Alpha value
ALPHA=0.15

## set line color
linecol="green"

####
####
####
####

## get expected values for DGRP data
null.l=fread("null_subset.pnps",header=F)
colnames(null.l)<-c("x","Chrom","NS","SS","y")
NL<-filter(null.l,Chrom == "genomewide")

NL.0 <- filter(NL,x==0)$y

## at first read SNAPE pnps data and filter genomewide values
DATA.pools=fread("PoolSNP_full.pnps",
                      header=T,
                      stringsAsFactors = F)

###
###
###

DATA.pools.group <- DATA.pools%>%
  filter(Chrom =="genomewide") %>%
  group_by(MAC,Chrom) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )

####

## keep pnps without MAF filtering only
DATA.pools.group.MAC50 <- DATA.pools %>%
  filter(MAC == 50 & Chrom =="genomewide")
#DATA.snape.group.MAF0

#read CSV of private SNPs
DATA.ps=read.csv("PoolSNP_full.ps",
                 header=T,
                 stringsAsFactors = F,
                 sep = "\t")

DATA=merge(DATA.ps,DATA.pools.group.MAC50, by.x="POP",by.y="POP")

###

DATA$private=log10(DATA$N)

## calculated Mean/SD and threshold based on Mean+2SD for pNpS data
Mean.pNpS=mean(DATA$pNpS)
SD.pNpS=sd(DATA$pNpS)
th.pNpS=Mean.pNpS+1.96*SD.pNpS

###

DATA$TH.pNpS <-DATA$pNpS
DATA$TH.pNpS[DATA$TH.pNpS<th.pNpS]<-NA
DATA$TH.pNpS[DATA$TH.pNpS>=th.pNpS]<-"Exclude"
DATA$TH.pNpS[is.na(DATA$TH.pNpS)]<-"Keep"

#save(DATA, file = "PNSPS.Keep.remove.Rdata")
DATA <- get(load("PNSPS.Keep.remove.Rdata"))
exclude.keep.pnps = DATA[,c(1,7,8,9)]
names(exclude.keep.pnps)[1] = "sampleId"

###
full_join(exclude.keep.pnps,
qc.dat.miss.cov.etc) %>% 
  mutate(TH.Cov = case_when(Cov > 10 ~ "Keep",
                            Cov <= 10 ~ "Exclude")) %>% 
  mutate(TH.Miss = case_when(Miss > 0.3 ~ "Exclude",
                              Miss <= 0.3 ~ "Keep")) %>% 
  mutate(TH.SimCont = case_when(SimCont > 0.15 ~ "Exclude",
                                SimCont <= 0.15 ~ "Keep")) -> final.calls
    
  save(final.calls, file = "final.calls.QC.Rdata")



