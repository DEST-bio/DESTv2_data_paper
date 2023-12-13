# ijob -A berglandlab -c1 -p standard --mem=10G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(poolfstat)
  library(SeqArray)
  library(data.table)
  library(gtools)

### load in SNPs to use
  load("/standard/vol186/bergland-lab/DEST_v2/snps_to_use.Rdata")
  load(file="/standard/vol186/bergland-lab/DEST_v2/GLMER_output_EuropeSeasonality.Rdata")
  table(snps$variant.id==sort(snps$variant.id))
  setkey(snps, variant.id)
  table(snps$variant.id==sort(snps$variant.id))

### sample metadata
  samps <- fread("/standard/vol186/bergland-lab/DEST_v2/seasonal_sets_DESTv2.csv")
  samps2 <- fread("/standard/vol186/bergland-lab/Gio/dest_v2.samps_8Jun2023.csv")
  samps <- merge(samps, samps2, by="sampleId")

### genofile
  genofile <- seqOpen("/standard/vol186/bergland-lab/DEST_v2/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds", allow.duplicate=T)
  seqSetFilter(genofile, sort(snps$variant.id))   #### sort error file

### get allele frequencies
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

### subset to samples
  ad.eu <- ad[samps$sampleId, ]
  dp.eu <- dp[samps$sampleId, ]

  ad.eu[is.na(ad.eu)] <- 0
  dp.eu[is.na(dp.eu)] <- 0


### make  poolfstat object
  snp.meta <- data.frame(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        ref=seqGetData(genofile, "$ref"),
                        alt=seqGetData(genofile, "$alt"))

### slice poolfstat object
  foreach(i=1:50)%dopar%{
    #i<-1
    message(i)
    index <- seq(i, dim(snps)[1], by = 50)
    #subpool <- pooldata.subset(pool, snp.index = index, return.snp.idx=T)
    subpool <- new("pooldata",
                npools=dim(ad.eu[,index])[1], #### Rows = Number of pools
                nsnp=dim(ad.eu[,index])[2], ### Columns = Number of SNPs
                refallele.readcount=t(ad.eu[,index]),
                readcoverage=t(dp.eu[,index]),
                poolsizes=samps$nFlies * 2,
                poolnames = rownames(ad.eu),
                snp.info = snp.meta[index,])

    pooldata2genobaypass(subpool, writing.dir = "/standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/", prefix = paste0("subpool_", i))
  }


### make contrast object
  cont <- data.table(sampleId=seqGetData(genofile, "sample.id"))
  cont <- merge(cont, samps[,c("sampleId", "season"), with=F], by="sampleId")

  ### load in old contrast file
  cont[,old:=t(fread("/standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt"))[,1]]
  table(cont$season, cont$old)
  cont[,POP:=1:138]
  save(cont, file="/standard/vol186/bergland-lab/alan/dest_baypass/contrast_pop.Rdata")
