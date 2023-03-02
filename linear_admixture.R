# ijob -c 4 --mem=8G -p dev -A berglandlab_standard

# module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)

### open GDS
  setwd("/project/berglandlab/DEST/gds")
  genofile <- seqOpen("dest.all.PoolSNP.001.5.25Feb2023.norep.ann.gds")
  length(seqGetData(genofile, "sample.id"))
  length(seqGetData(genofile, "variant.id"))

### open metadata
  setwd("/scratch/aob2x/DESTv2/")
  samps <- fread("populationInfo/dest_v2.samps_25Feb2023.csv")

### define test pops
  af_pop <- samps[continent=="Africa"][nFlies>100]$sampleId
  eu_pop <- sample(samps[continent=="Europe"]$sampleId, 1)

  test_pop <- samps[continent=="North_America"][grepl("Nor", locality)][1]$sampleId


### snp library
  snps.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        chr=seqGetData(genofile, "chromosome"))


### get allele frequencies
  ### define test pops
    af_pop <- samps[continent=="Africa"][nFlies>100]$sampleId
    eu_pop <- sample(samps[continent=="Europe"]$sampleId, 1)

    ### NC
      test_pop <- samps[continent=="North_America"][grepl("Nor", locality)][1]$sampleId

  ### seq set filter
    seqSetFilter(genofile, variant.id=sample(snps.dt[nAlleles==2]$variant.id, 100), sample.id=c(af_pop, eu_pop, test_pop))

  ### build model
    dat <- t(seqGetData(genofile, "annotation/format/FREQ")$data)
    dat <- as.data.table(dat)
    setnames(dat, names(dat), seqGetData(genofile, "sample.id"))
    setnames(dat, names(dat), c("eu", "us", "af"))

    summary(lm(us~0+eu+af, dat))


    ### AUS
      test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
    ### seq set filter
      seqSetFilter(genofile, variant.id=sample(snps.dt[nAlleles==2]$variant.id, 1000), sample.id=c(af_pop, eu_pop, test_pop))

    ### build model
      dat<- t(seqGetData(genofile, "annotation/format/FREQ")$data)
      dat <- as.data.table(dat)
      setnames(dat, names(dat), seqGetData(genofile, "sample.id"))

      setnames(dat, names(dat), c("aus", "eu", "af"))

      summary(lm(aus~0+eu+af, dat))










  dat <- expand.grid(dat.w)[,1]

  dat <- data.table(af_pop=seqGetData(genofile, "annotation/format/FREQ")$data,
                    eu_pop=seqGetData(genofile, "annotation/format/FREQ")$data,
                    test_pop=seqGetData(genofile, "annotation/format/FREQ")$data)
