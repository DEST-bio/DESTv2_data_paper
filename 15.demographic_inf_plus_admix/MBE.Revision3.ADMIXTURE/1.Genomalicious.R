#### Moments admixture
# module load Rgeospatial
library(genomalicious)
library(SeqArray)
library(foreach)
library(sp)
##pairsF="/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/pairs_guide_file.txt"

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
pairsF=args[2]

pairs = fread(pairsF)
#1.Create output the directory
system(paste("mkdir", "dadi_objects", sep = " "))
system(paste("mkdir", "L_meta_objects", sep = " "))

#2. Load the pairs data
SAMPLE_INPUT=c(pairs$Parent1[i],pairs$Parent2[i],pairs$Derived[i])

#### Read in GDS file
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")
meta <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")
# make SNP table - old
seqResetFilter(genofile)
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
### convert to format for genomalicious
tf <- TRUE
dat <- as.data.table(dat)
dat[,locus:=seqGetData(genofile, "variant.id")[tf]]
dat[,ref:=seqGetData(genofile, "$ref")[tf]]
dat[,alt:=seqGetData(genofile, "$alt")[tf]]

datl <- melt(dat, id.vars=c("locus", "ref", "alt"))
setnames(datl, c("variable", "value"), c("POOL", "FREQ"))
######

### Calculate NEff
t(dp) %>%
  data.frame() -> dp.df
names(dp.df) = seqGetData(genofile, "sample.id")
colMeans(dp.df, na.rm = T) %>%
  data.frame(DP=.) %>%
  mutate( sampleId=rownames(.)) ->
  dp.df
####


####
meta %>%
  filter(sampleId %in% seqGetData(genofile, "sample.id")) %>%
  select(sampleId, nFlies) -> inds_nfly

left_join(dp.df, inds_nfly) %>%
  mutate(ne = floor((nFlies*DP)/(nFlies+DP-1)) )->neff

datl <- merge(datl, neff[,c("sampleId", "ne")], by.x="POOL", by.y="sampleId")
datl[,FREQ:=1-FREQ]
datl <- na.omit(datl)

datl %>%
  filter(POOL == pairs$Parent1[i]) %>%
  mutate(POOL = "A_parent1") ->
  datl1

datl %>%
  filter(POOL == pairs$Parent2[i]) %>%
  mutate(POOL = "B_parent2") ->
  datl2

datl %>%
  filter(POOL == pairs$Derived[i]) %>%
  mutate(POOL = "C_derived") ->
  datl3


rbind(datl1, datl2, datl3) -> datl.relab

sfs_method="probs"

dadi <- dadi_inputs(
  datl.relab,
  type = "freqs",
  popCol = "POOL",
  locusCol = "locus",
  refCol = "ref",
  altCol = "alt",
  freqCol = "FREQ",
  indsCol = "ne",
  freqMethod = sfs_method
)
dadi <- na.omit(dadi) 

L=dim(dadi)[1]
####
### Write dadi file

fn <- paste("dadi_objects/",
            sfs_method, ".",
            pairs[i,]$Parent1, ".",
            pairs[i,]$Parent2, ".",
            pairs[i,]$Derived, ".",
            "delim", sep="")

write.table(dadi,
            file=fn,
            sep="\t", quote=F, row.names=F)

### Metadata
neff %>%
  filter(sampleId == pairs[i,]$Parent1) %>% .$ne -> an1_ne
neff %>%
  filter(sampleId == pairs[i,]$Parent2) %>% .$ne -> an2_ne
neff %>%
  filter(sampleId == pairs[i,]$Derived) %>% .$ne -> Der_ne

fn2 <- paste("L_meta_objects/",
            sfs_method, ".",
            pairs[i,]$Parent1, ".",
            pairs[i,]$Parent2, ".",
            pairs[i,]$Derived, ".",
            "meta", sep="")


write.table(data.frame(an1_ne, an2_ne, Der_ne, L),
            file=fn2,
            sep="\t", quote=F, row.names=F)

