#### DEST Analysis FST across time and across space
#### Script 1 -- find the samples used for the FST analysis.

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
SEL= as.numeric(args[1]) 
##1-8

root="/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/outfiles/"

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(foreach)
library(gtools)

library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

ingds = "/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"

##### step 1 
### filter samps by pass
#### Filter by cluster
#### Cluster 2 == EUW
#### Cluster 8 == EUE

samps %>%
filter(Recommendation == "Pass")  ->
samps.pass

#### First calculate temporal FST
#samps.pass %>%
#filter(year >= 2016) %>%
#group_by(locality, year) %>%
#summarize(N=n()) %>%
#reshape2::dcast(locality~year, value.var = "N", fill= 0) %>%
#mutate(tot=`2016`+`2017`+`2018`+`2019`+`2020`+`2021`
#       ) %>%
#filter(tot>=5) %>%
#.$locality -> select_locs

select_locs = c("CH_Vau_Vul",
"DE_Bay_Mun",
"DK_Mid_Kar",
"FI_Pir_Aka",
"TR_Ank_Yes",
"UA_Kie_Vys",
"UA_Ode_Ode",
"US_Vir_Cha")

#####
##spatial_samps_pery =
##foreach(i = select_locs,
##.combine="rbind"
##)%do%{
##
##samps.pass %>%
##filter(year > 2012) %>%
##group_by(year) %>%
##summarize(N=n()) -> 
##n.y.samps
##
##
##for(k in 1:dim(n.y.samps)[1]){
##y.n.max = n.y.samps$year[which(n.y.samps$N == 
##sort(n.y.samps$N, decreasing = T)[k]) ]
##
##if(i %in% samps.pass[year == y.n.max]$locality){
##message("yes")
##o.l = (data.frame(i,k, y.n.max ))
##break
##}else{
##message("no")
##} ### if
##} ### for loop in
##
##return(o.l)
##
##} ### foreach loop out
##
#####
#####
#####
#####
##### Begin temporal code
samps.pass %>%
filter(locality == select_locs[SEL]) %>%
mutate(collectionDate = paste(
min_day,
min_month,
year,
sep = "-"
)) %>% mutate(collectionDate = as.Date(collectionDate, 
format = "%d-%m-%Y")) ->
working.obj

L = dim(working.obj)[1]


comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

print("Create combination vector")

comp_vector %<>%
  as.data.frame() %>%
  mutate(day_diff = NA)


for(i in 1:dim(comp_vector)[1]) {
  
  date1=working.obj$collectionDate[comp_vector[i,1]]
  date2=working.obj$collectionDate[comp_vector[i,2]]
  
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  
  comp_vector$pop1[i] = working.obj$city[comp_vector[i,1]]
  comp_vector$pop2[i] = working.obj$city[comp_vector[i,2]]
  
  comp_vector$samp1[i] = working.obj$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = working.obj$sampleId[comp_vector[i,2]]
  
}

#####
########################
### open GDS file
genofile <- seqOpen(ingds)

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=working.obj$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

snps.dt %<>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))
seqSetFilter(genofile, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
### select sites
  seqSetFilter(genofile,
              snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05][af>.2]$variant.id)
              
              
### get allele frequency data
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

###################
### Next part
### Calculate FST

#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector)[1]),
  samp2 = rep(NA, dim(comp_vector)[1]),
  FST = rep(NA, dim(comp_vector)[1])
)

###
###
###
for(i in 1:dim(comp_vector)[1]){
  
  print(i/dim(comp_vector)[1] * 100)
  
  samps_to_compare = c(comp_vector$samp1[i], comp_vector$samp2[i])
  
  pool_sizes = c(working.obj$nFlies[which(working.obj$sampleId == comp_vector$samp1[i])],
                 working.obj$nFlies[which(working.obj$sampleId == comp_vector$samp2[i])])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)
  
  
  fst.out <- computeFST(pool, method = "Anova")
  
  outfile$samp1[i] = comp_vector$samp1[i]
  outfile$samp2[i] = comp_vector$samp2[i]
  outfile$FST[i] = fst.out$FST
  
}## i   
            
            
left_join(comp_vector, outfile) -> Out_comp_vector_samepops

file.name = paste(select_locs[SEL],"time",sep = ".")

save(Out_comp_vector_samepops,
file = paste(root,file.name, ".Rdata", sep = ""))


