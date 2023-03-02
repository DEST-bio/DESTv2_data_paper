# ijob -c 4 --mem=8G -p dev -A berglandlab_standard

# module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

### libraries
  library(SeqArray)
  library(data.table)

### open GDS
  setwd("/project/berglandlab/DEST/gds")
  genofile <- seqOpen("dest.all.PoolSNP.001.5.25Feb2023.norep.ann.gds")
  length(seqGetData(genofile, "sample.id"))
  length(seqGetData(genofile, "variant.id"))


  genofile2 <- seqOpen("dest.all.PoolSNP.001.50.10Nov2020.ann.gds")
  length(seqGetData(genofile2, "sample.id"))
  length(seqGetData(genofile2, "variant.id"))

### open metadata
  setwd("/scratch/aob2x/DESTv2/")
  samps <- fread("populationInfo/dest_v2.samps_25Feb2023.csv")

### depth
  snps.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        chr=seqGetData(genofile, "chromosome"))

  seqSetFilter(genofile, variant.id=sample(snps.dt[nAlleles==2]$variant.id, 1e5))


  ac <- seqGetData(genofile, "annotation/info/AC")









  rd.dt <-

  dt.ag <- dt[,list(nMissing=sum(is.na(x)), mrd=(mean(x, na.rm=T))), sampleId]

  dt.ag <- merge(dt.ag, samps, by="sampleId")
  dt.ag[,nEff:=(mrd*2*nFlies)/(mrd+2*nFlies)]
  summary(lm(mrd~set, dt.ag))

  dt.ag[,list(mu=median(nEff, na.rm=T)), list(set)]
  dt.ag[,list(mu=median(nEff, na.rm=T)), list(continent)]


### who is missing
  setkey(samps, sampleId)
  samps2 <- merge(samps, data.table(sampleId=seqGetData(genofile, "sample.id"), mapped=T, key="sampleId"), all.x=T, all.y=T)
  samps2[is.na(mapped), mapped:=F]

  table(samps2$mapped, samps2$set)
