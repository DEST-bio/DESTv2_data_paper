library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(foreach)
library(gtools)
library(doParallel)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(doMC)
registerDoMC(4)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")
ingds = "/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"

####
load("y1comps.guide.Rdata")
setDT(y1comps)

### Divide into 32 jobs
jobIds=seq(from = 1, to =29513, by = 32 )

root="/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/outfiles_1y/"
args = commandArgs(trailingOnly=TRUE)
SEL= as.numeric(args[1]) 
start.i = jobIds[SEL]

##1-923



#####
y1comps[c(start.i:(start.i+31)),] ->
  tmp.select
  
#samples.select = unique(unique(tmp.select$samp1),unique(tmp.select$samp2))
##### Begin temporal code
#####
########################
### open GDS file
genofile <- seqOpen(ingds)

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

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
  samp1 = rep(NA, dim(tmp.select)[1]),
  samp2 = rep(NA, dim(tmp.select)[1]),
  FST = rep(NA, dim(tmp.select)[1])
)

###
###
###
for(i in 1:dim(tmp.select)[1]){
  
  print(i/dim(tmp.select)[1] * 100)
  
  samps_to_compare = c(tmp.select$samp1[i], tmp.select$samp2[i])
  
  pool_sizes = c(
    samps$nFlies[which(samps$sampleId == tmp.select$samp1[i])],
    samps$nFlies[which(samps$sampleId == tmp.select$samp2[i])]
    )
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)
  
  
  fst.out <- computeFST(pool, method = "Anova")
  
  outfile$samp1[i] = tmp.select$samp1[i]
  outfile$samp2[i] = tmp.select$samp2[i]
  outfile$FST[i] = fst.out$FST
  
}## i   


left_join(tmp.select, outfile) -> Out_comp_vector_samepops

file.name = paste("job", jobIds[SEL],"y1comp",sep = ".")

save(Out_comp_vector_samepops,
     file = paste(root,file.name, ".Rdata", sep = ""))
