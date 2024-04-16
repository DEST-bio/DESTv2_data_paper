### Evaluate Italy
library(tidyverse)
library(data.table)
library(SeqArray)
library(magrittr)

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
seqSetFilter(genofile, 
             sample.id=filter(samps, sampleId %in% samps.for.model.train)$sampleId,    
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
model2.0$DAPC$ind.coord %>% rownames -> samps.for.model.train


New_samples=dat_DEST1_markers.naImp %>%
  .[which(rownames(.) %in% samps.dest.samps.NewSamps$sampleId),]

predict.dapc(model2.0$DAPC, newdata=New_samples) -> DEST1_predictions
