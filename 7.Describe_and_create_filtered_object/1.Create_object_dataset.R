### DESCRIBE DATASET
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)
library(matrixStats)
library(foreach)
library(doParallel)
library(doMC)
registerDoMC(4)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)

library(DescTools)

###
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))

#####
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- vroom("dest_v2.samps_25Feb2023.csv")

samps %>% 
  filter(set != "dgn") ->
  samps.flt

####

filtering.dt %>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  chosen.snps
####

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt

####
snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05]$variant.id)

#### extract matrices of cov and calls
### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

### divide the matrices to generate allele frequency calls
dat <- ad/dp

#####
#####
### check the dimensions of the allele frequency matrix
dim(dat)  

####
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

rownames(dat) <- seqGetData(genofile, "sample.id")
rownames(ad) <- seqGetData(genofile, "sample.id")
rownames(dp) <- seqGetData(genofile, "sample.id")

####
####
####

dat %>%
  t() %>% 
  as.data.frame -> dat_t

dat_t[,samps.flt$sampleId] -> dat_t.flt

dat_t.flt %>%
  mutate(SNP_id.SNP_num = rownames(.)) %>% 
  separate(SNP_id.SNP_num, into = c("chr", "pos", "SNP_id")) ->
  dat_t.flt.met

#####

###
## Some characterizations of AFs and subsequent filtering
MeanAF=c()
MinAF=c()

apply(dat_t.flt.met[,-which(colnames(dat_t.flt.met) %in% c("chr", "pos", "SNP_id"))],
      1, FUN=mean, na.rm=TRUE ) -> MeanAF
data.frame(SNP_id = dat_t.flt.met$SNP_id, MeanAF) -> MeanAF

apply(dat_t.flt.met[,-which(colnames(dat_t.flt.met) %in% c("chr", "pos", "SNP_id"))],
      1, FUN=min, na.rm=TRUE ) -> MinAF
data.frame(SNP_id = dat_t.flt.met$SNP_id, MinAF) -> MinAF

cbind(dat_t.flt.met, MeanAF[-1], MinAF[-1]) -> dat_t.met.ag

##
dat_t.met.ag %>%
  #.[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MeanAF > 0.01),] ->  ### This samples only polymorphic sites
  dat_t.met.ag.Mean0.01

#####
#####
#####

count_NA = function(x){
  return(sum(is.na(x)))
}
MissDat=c()

pool_cols = grep("SNP_id|MeanAF|MinAF", colnames(dat_t.met.ag.Mean0.01), invert = T)

apply(dat_t.met.ag.Mean0.01[,pool_cols],
      1, FUN=count_NA ) -> MissDat

n_pools = length(pool_cols)

data.frame(missing_rate = c(MissDat/n_pools) ) -> MissDat

cbind(dat_t.met.ag.Mean0.01, MissDat) -> dat_t.met.ag.Mean0.01.miss

dat_t.met.ag.Mean0.01.miss %>%
  as.data.frame() %>%
  filter(missing_rate < 0.01) ->  
  dat_t.met.ag.Mean0.01.miss0.01

#####
dat_t.met.ag.Mean0.01.miss0.01 %>% 
  group_by(chr) %>%
  summarize(N=n())

####
dat_t.met.ag.Mean0.01.miss0.01 %>% 
  rownames() -> flt.snps.vector

####
dat[,flt.snps.vector] ->
  dat.o

save(dat.o, file = "DEST2.flt.meanAF0.01.miss0.01.Rdata" )
save(dat.o, file = "/project/berglandlab/DEST/SNP_Filtering_files/DEST2.flt.meanAF0.01.miss0.01.Rdata" )





