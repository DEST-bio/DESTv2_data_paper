##### Calculate F3 statistics in pool-seq
##### 


args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1]) ## 1-347


### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
library(vroom)
library(poolfstat)


### open GDS

genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.5.25Feb2023.norep.ann.gds")
length(seqGetData(genofile, "sample.id"))
length(seqGetData(genofile, "variant.id"))

### open metadata
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")
samps <- fread("dest_v2.samps_25Feb2023.csv")

### snp library
snps.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqGetData(genofile, "$num_allele"),
                      chr=seqGetData(genofile, "chromosome"))

### define test pops
af_pop <- samps[continent=="Africa"][nFlies>100]$sampleId
eu_pop <- samps[continent=="Europe"]$sampleId[k]
AUS_samps = samps[continent=="Oceania"]$sampleId

########
au_pop = AUS_samps[i]

seqSetFilter(genofile, 
             variant.id=sample(snps.dt[nAlleles==2]$variant.id, 10000),
             sample.id=c(af_pop,eu_pop,au_pop))

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")



####
poolsizes = c(
  filter(samps, sampleId == rownames(dp)[1])$nFlies,
  filter(samps, sampleId == rownames(dp)[2])$nFlies,
  filter(samps, sampleId == rownames(dp)[3])$nFlies
)

snp.info = data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position")
                      )

pool <- new("pooldata",
            npools=dim(ad)[1], #### Rows = Number of pools
            nsnp=dim(ad)[2], ### Columns = Number of SNPs
            refallele.readcount=t(ad),
            readcoverage=t(dp),
            poolsizes=poolsizes * 2,
            poolnames = rownames(dp),
            snp.info = snp.info
            )
## snp.info A data frame (nsnp rows and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, Reference allele name and Alternate allele name (if available)
##
fstats.dat <- 
compute.fstats(
  pool,
  nsnp.per.bjack.block = 100,
  computeDstat = TRUE,
  verbose = TRUE
)

fstats.dat@f3star