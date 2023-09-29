### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)
library(data.table)
library(FactoMineR)
library(adegenet)
library(zoo)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)


setwd("/scratch/yey2sn/DEST2_analysis/dapc_dims")

gims.snps <- get(load("DEST.2.0.GIMS.Rdata"))
#####
message("now loading AF")
####

meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"

samps <- fread(meta_git)
setDT(samps)

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps.o.nsim.noBeij


###
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)

seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

snps.dt %>%
filter(SNP_id %in% gims.snps$SNP_id ) -> 
snps.dt.gims

###
### Remove singletons
samps.o.nsim.noBeij %>% 
  group_by(province) %>%
  summarize(N = n()) %>%
  filter(N == 1) -> singleton.cities

samps.o.nsim.noBeij %>%
  filter(!province %in% singleton.cities$province) ->
  gim.samps

gim.samps$province = gsub("Kiev City", "Kiev", gim.samps$province)

####
  seqSetFilter(genofile, 
               sample.id=gim.samps$sampleId, 
               variant.id=snps.dt.gims$variant.id)
####
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  ad.matrix = ad$data
  dat <- ad.matrix/dp
  dim(dat)
  colnames(dat) <- paste(seqGetData(genofile, "chromosome"),
                         seqGetData(genofile, "position") 
                         , sep="_")
  rownames(dat) <- seqGetData(genofile, "sample.id")

### clean AF tables from NA
#apply(gims.flt.AF, 2, function(x) sum(is.na(x)) ) -> Count.of.NAs
message("now na.aggregate")

dat_naImp = na.aggregate(dat)


###
message("now training xval")
xval_province <- xvalDapc(dat_naImp, 
                          as.factor(gim.samps$province), 
                          n.pca.max = 300, 
                          training.set = 0.9,
                          result = "groupMean", 
                          center = TRUE, 
                          scale = FALSE, 
                          n.pca = NULL, 
                          n.rep = 30, 
                          xval.plot = FALSE)

save(xval_province, file = "xval_province.DEST2.0.Rdata")

