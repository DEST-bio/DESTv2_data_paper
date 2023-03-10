#### Merging Protocol for Australian flies

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

#####
##########


system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- fread("dest_v2.samps_25Feb2023.csv")

##samps %>% filter(country == "Germany") %>% head

samps %>% filter(!is.na(loc_rep) & collector == "Fournier-Level et al") ->
  samps_to_merge

#########
### Load the genotype locality
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds", allow.duplicate=T)
seqResetFilter(genofile)
seqSetFilter(genofile, 
             sample.id=samps_to_merge$sampleId)
#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/project/berglandlab/DEST/SNP_Filtering_files/DESTv2.SNPmeta.filter.Rdata"))
filtering.dt %<>%
  filter(is.na(libs))
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))
snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]
####
### Begin Merging

samps_to_merge$loc_rep %>% unique -> iths

Merge.ADs =
foreach(i = iths, .combine = "rbind")%do%{
  seqResetFilter(genofile)
  seqSetFilter(genofile, 
               sample.id=filter(samps_to_merge, loc_rep == i )$sampleId, 
               variant.id=filtering.dt$id)


  ad <- seqGetData(genofile, "annotation/format/AD")$data
  colnames(ad) <- paste(seqGetData(genofile, "chromosome"), 
                         seqGetData(genofile, "position") , 
                         paste("snp", seqGetData(genofile, "variant.id"), sep = ""), 
                         sep="_")
  
  ad.sum = colSums(ad, na.rm = T)
  data.frame(t(as.data.frame(ad.sum))) -> ad.sum.dft
  rownames(ad.sum.dft) = gsub(" ", "_", i)
  colnames(ad.sum.dft) = gsub("^X", "", colnames(ad.sum.dft))
  
  return(ad.sum.dft)
  
}

########
########
########
########
Merge.DPs =
  foreach(i = iths, .combine = "rbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=filter(samps_to_merge, loc_rep == i )$sampleId, 
                 variant.id=filtering.dt$id)
    
    
    dp <- seqGetData(genofile, "annotation/format/DP")
    colnames(dp) <- paste(seqGetData(genofile, "chromosome"), 
                          seqGetData(genofile, "position") , 
                          paste("snp", seqGetData(genofile, "variant.id"), sep = ""), 
                          sep="_")
    
    dp.sum = colSums(dp, na.rm = T)
    data.frame(t(as.data.frame(dp.sum))) -> dp.sum.dft
    rownames(dp.sum.dft) = gsub(" ", "_", i)
    colnames(dp.sum.dft) = gsub("^X", "", colnames(dp.sum.dft))
    
    return(dp.sum.dft)
    
  }

## --> make dt
samps_to_merge %>%
  group_by(loc_rep) %>%
  slice_head() %>%
  mutate(sampleId = gsub(" ", "_", loc_rep)) ->
  merged_metadata
#####


#### ---> Regular DP/AD
#### ---> Regular DP/AD
#### ---> Regular DP/AD
#### ---> Regular DP/AD

samps %>% filter(!(sampleId %in% samps_to_merge$sampleId)) ->
  regular.samps

dim(regular.samps)[1] + dim(samps_to_merge)[1] ==  dim(samps)[1]

#########
### Load the genotype locality
seqResetFilter(genofile)
seqSetFilter(genofile, 
             sample.id=regular.samps$sampleId,
             variant.id=snps.dt$variant.id)
ad <- seqGetData(genofile, "annotation/format/AD")$data
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), 
                      seqGetData(genofile, "position") , 
                      paste("snp", seqGetData(genofile, "variant.id"), sep = ""), 
                      sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")
dp <- seqGetData(genofile, "annotation/format/DP")
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), 
                      seqGetData(genofile, "position") , 
                      paste("snp", seqGetData(genofile, "variant.id"), sep = ""), 
                      sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")
####
#### Sanity Checks!
dim(ad)[2] + dim(Merge.ADs)[2] ==  dim(dp)[1] + dim(Merge.DPs)[1]

intersect(colnames(ad),colnames(Merge.ADs)) -> shared_SNPs

joint.ad = rbind(as.data.frame(ad)[shared_SNPs], Merge.ADs[shared_SNPs])
joint.dp = rbind(as.data.frame(dp)[shared_SNPs], Merge.DPs[shared_SNPs])

##### Create metadata
rbind(samps, merged_metadata) -> joint.metadata

### SAVE OBJECTS
### 
save(joint.ad,
     file = "/project/berglandlab/DEST2.0_working_data/joint.ad.matrix.Rdata")
save(joint.dp,
     file = "/project/berglandlab/DEST2.0_working_data/joint.dp.matrix.Rdata")
save(joint.metadata,
     file = "/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata")

#### CREATE AF object.
#### 
joint.dat = joint.ad/joint.dp
save(joint.dat,
     file = "/project/berglandlab/DEST2.0_working_data/joint.AFs.Rdata")

###


