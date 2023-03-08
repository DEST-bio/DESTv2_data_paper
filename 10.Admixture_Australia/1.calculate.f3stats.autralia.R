##### Calculate F3 statistics in pool-seq
##### 



library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(zoo)
library(adegenet)
library(reshape2)
library("geosphere")
library(gmodels)
library(ggforce)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

##3
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))

###
###
#####
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

samps <- vroom("dest_v2.samps_25Feb2023.csv")
grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos
samps[-c(sim.pos, Beijing.pos),] -> samps

samps %>%
  filter(continent %in% c("Europe","Africa","Oceania") ) ->
  eu.afr.oce.samps

####
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=eu.afr.oce.samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

###
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
dat <- ad$data/dp
dim(dat)



