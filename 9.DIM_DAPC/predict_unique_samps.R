### Evaluate Italy
library(tidyverse)
library(data.table)
library(SeqArray)
library(magrittr)
library(adegenet)

gims.snps <- get(load("DEST.2.0.GIMS.Rdata"))
###
genofile <- seqOpen("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds", allow.duplicate=T)
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

snps.dt %>%
  filter(SNP_id %in% gims.snps$SNP_id ) -> 
  snps.dt.gims

####
samp_interest = "IT_Sas_Rec_1_2018-10-18"
  
seqSetFilter(genofile, 
             sample.id=filter(samps, sampleId %in% samp_interest)$sampleId,    
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


####
model2.0 <- get(load("xval_province.DEST2.0.Rdata"))

predict.dapc(model2.0$DAPC, newdata=dat) -> DEST2_predictions

t(DEST2_predictions$posterior) %>% data.frame() -> o
o %>% arrange(-IT_Sas_Rec_1_2018.10.18)
