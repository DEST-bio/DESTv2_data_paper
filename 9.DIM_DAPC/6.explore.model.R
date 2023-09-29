### Explore Model of DEST 2.0
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

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(forcats)
library(FactoMineR)

####
####
####
#sftp://rivanna.hpc.virginia.edu/scratch/yey2sn/DEST2_analysis/dapc_dims/xval_province.DEST2.0.Rdata

gims.snps <- get(load("DEST.2.0.GIMS.Rdata"))

model2.0 <- get(load("xval_province.DEST2.0.Rdata"))
model2.0$DAPC$ind.coord %>% rownames -> samps.for.model.train

#model2.0$DAPC$pca.loadings %>% rownames() %>% sort -> GIMs

#####
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv"

samps <- fread(meta_git)
setDT(samps)

#samps %>% filter(sampleId %in% samps.for.model.train) ->
#samps.used.in.modeltrain 
####

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
samps.for.model.train

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



#### ---> Creating global object 
dat_naImp = na.aggregate(dat)

#### Part 1. Split sample set..
#### 
#i=1
out2 = list()
out =
foreach(i=1:dim(dat_naImp)[1],
        #i=1:10,
        .combine = "rbind",
        .errorhandling = "remove"
        )%do%{
          
          message(i)

dat_naImp[-which(rownames(dat_naImp) == samps$sampleId[i]),] -> dat.gims_naImp_loo

dat_naImp[samps$sampleId[i],] %>%
  as.data.frame %>%
  t() ->
  dat.gims_naImp_anchor

###
samp.train =
  filter(samps, sampleId %in% samps.for.model.train) %>% 
  filter(!sampleId == samps$sampleId[i])

##########
##########
dapc(dat.gims_naImp_loo, 
     grp=samp.train$province, 
     n.pca=100, n.da=75) ->
  model.lo_out

predict.dapc(model.lo_out, 
newdata=dat.gims_naImp_anchor) -> LOO_predictions

  data.frame(pop = rownames(t(LOO_predictions$posterior)),
         P = t(LOO_predictions$posterior)) ->
    tmp.p
names(tmp.p)[2] = "P"
  
tmp.p %>%
  slice_max(P) ->
  post.loo

predicted.pop = post.loo$pop
predicted.post = post.loo$P

samps[sampleId == samps$sampleId[i]] %>%
  group_by(province) %>%
  summarize(real.lat = mean(lat),
            real.long = mean(long),
            cont = unique(continent)
            ) -> real.mean

samps[province == post.loo$pop] %>%
  group_by(province) %>%
  summarize(pred.lat = mean(lat),
            pred.long = mean(long),
  ) -> pred.mean

cbind(predicted.pop, 
      predicted.post, 
      pred.mean$pred.lat, 
      pred.mean$pred.long, 
      real.mean, sampleid = samps$sampleId[i]) -> o
        
      o %>% 
  group_by(sampleid) %>%
  mutate(hav_d = distHaversine(
    matrix(c(`pred.mean$pred.long`, real.long, `pred.mean$pred.lat`, real.lat),nrow = 2))/1000) ->
  o
   
out2[[i]] = o
return(o)

        }

do.call("rbind", out2) -> gims.out.v2

save(gims.out.v2, file = "DEST2.0.model.predictions.Rdata")

###
