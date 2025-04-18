#### Moments admixture
# module load Rgeospatial
library(genomalicious)
library(SeqArray)
library(foreach)
library(sp)

####
args = commandArgs(trailingOnly=TRUE)

job=as.numeric(args[1])
popset=args[2]
dir= args[3]
Meta_dir=args[4]
DEST_gds_dir=args[5]

print(c(job,popset,dir,Meta_dir,DEST_gds_dir))

SAMPLE_INPUT=c("US_Tex_Sto_10_2013-09-15", 
               "US_Tex_Sto_1_2013-09-15", 
               "US_Tex_Sto_2_2013-09-15")




#### Read in GDS file
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")

# make SNP table - old
message("making snp table")
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))
## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

message("get poly")
setkey(snps.dt, chr)
seqSetFilter(genofile, 
             sample.id=SAMPLE_INPUT,
             variant.id=snps.dt[J(c("2L", "2R", "3L", "3R"))]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
dat <- ad$data/dp

### parse
dim(dat)
rownames(dat) <- seqGetData(genofile, "sample.id")
dat <- t(dat)

####
tf<-TRUE
##  dim(dat)
##}


### Calculate Neff
### load average effective read depth
setkey(dep, sampleId)
setkey(samps, sampleId)
neff <- merge(dep[J(colnames(dat))], samps[J(colnames(dat))], by="sampleId")[,c("sampleId", "nFlies", "mu.25"), with=F]
neff[,ne:=round((2*nFlies*mu.25)/(2*nFlies+mu.25))]
neff[,ne:=floor(ne/2)] ### note, this is the number of diploid individuals...



