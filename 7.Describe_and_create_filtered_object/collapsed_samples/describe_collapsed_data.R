### Characterize Data set
### 

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
library(forcats)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)


#####
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.qc_merge.csv"

samps <- fread(meta_git)
setDT(samps)
#samps = samps[set!="dgn"]
samps = samps[!is.na(collapsedSamples)]

##### DESCRIBE ORIGINAL GDS
##### DESCRIBE ORIGINAL GDS

#### ----> Generate Basic Stats
####
####
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))
filtering.dt %<>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) 


#grep( "SIM" , samps$sampleId) -> sim.pos
#grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos
#samps[-c(sim.pos, Beijing.pos),] -> samps

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

snps.dt <- snps.dt[nAlleles==2]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %<>% filter(SNP_id %in% filtering.dt$SNP_id)
snps.dt %>% dim

####
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

########
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
ad.matrix = ad$data
dat <- ad.matrix/dp
dim(dat)

colnames(dp) <- paste(seqGetData(genofile, "chromosome"),
                      seqGetData(genofile, "position") 
                      , sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

#### Estimate NA %

rowSums(is.na(dp))/dim(dp)[2] -> Missing.data.calc
rowMeans(dp, na.rm = T) -> Means.Cov

####### COVERGAE DATA
Means.Cov.df =
  data.frame(Value = Means.Cov,
             Var = "Cov") %>% 
  mutate(sampleId = names(Means.Cov))# %>%
  #left_join(samps)

### Missing data
Missing.data.df =
  data.frame(Value = Missing.data.calc,
             Var = "Miss") %>% 
  mutate(sampleId = names(Missing.data.calc)) #%>%
  #left_join(samps)
####

#save(Miss.Cov.joint, file = "Miss.and.Cov.joint.Rdata")

### PCR DUPLICATES
### 
fns <- system("find /project/berglandlab/DEST/dest_mapped/DEST_plus_collapsed -name '*duplicates_report.txt'", intern=T)
pcr <- foreach(fn=fns, .combine = "rbind")%dopar%{
  #fn <- fns[1]
  message(fn)
  tmp <- fread(fn, skip="LIBRARY", nrows=1)
  data.table(pcrDup=tmp$PERCENT_DUPLICATION, READ_PAIRS_EXAMINED=tmp$READ_PAIRS_EXAMINED, sampleId=tstrsplit(fn, "/")[[7]])
}

pcr %>%
  as.data.frame() %>%
  dplyr::select(Value = pcrDup, sampleId) %>% 
  mutate(Value=Value, Var = "PCRdup")   ->
  pcr.mod

####

pcr.mod
Means.Cov.df
Missing.data.df

rbind(Means.Cov.df, Missing.data.df, pcr.mod) -> covs.miss.pcr

save(covs.miss.pcr, file = "/scratch/yey2sn/DEST2_analysis/filtering/collapsed_samps/collapsed.covs.miss.pcr.Rdata")
