# ijob -A berglandlab_standard -c5 -p standard --mem=40G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R
  args = commandArgs(trailingOnly=TRUE)
  jobId <- as.numeric(args[1]) # jobId=100

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(SeqArray)
  library(data.table)
  library(foreach)
  #library(tidyverse)
  #library(magrittr)
  #library(vroom)
  library(poolfstat)
  library(doMC)
  registerDoMC(5)

### metadata
  meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv"
  samps <- fread(meta_git)
  setnames(samps, "cluster2.0_k4", "cluster")
  samps[country=="Ukraine"]$cluster

### define tests
  samps <- samps[Recommendation=="Pass"]
  americas_aus <- unique(samps[grepl("America", continent) | grepl("Oceania", continent)]$locality)
  eurasia <- unique(samps[grepl("Europe", continent) | grepl("Asia", continent)]$locality)
  africa <- unique(samps[grepl("Africa", continent)]$locality)

  west_europe <- unique(samps[grepl("Europe", continent)][cluster==3]$locality)
  east_europe <- unique(samps[grepl("Europe", continent)][cluster==2]$locality)

  tests <- as.data.table(expand.grid(focal=americas_aus, eurasia, africa))
  tests2 <- as.data.table(expand.grid(focal=west_europe, east_europe, africa))

  tests <- rbind(tests, tests2)
  tests[,job:=rep(1:5000, each=ceiling(dim(tests)[1]/5000))[1:dim(tests)[1]]]
  tests <- tests[job==jobId]

### open GDS
  #genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)
  genofile <- seqOpen("/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")

### load SNP table
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=samps$sampleId)
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile),
                        ref_allele=seqGetData(genofile, "$ref"),
                        alt_allele=seqGetData(genofile, "$alt"))

  snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]
  snps.dt[,SNP_id:=paste(chr, pos, sep="_")]

#### add in annotation
  # snpdt.obj <- get(load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt.Rdata"))
  # setDT(snpdt.obj)
  # snpdt.obj %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
  # snpdt.obj.NoInv = snpdt.obj[invName == "none"]
  #
  # annotation <- get(load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_Cville_SNP_Metadata.Rdata"))
  # setDT(annotation)
  # names(annotation)[1:2] = c("chr","pos")
  # non.cod = c("intergenic_region","intron_variant","upstream_gene_variant","downstream_gene_variant")
  # annotation.non.cod = annotation[effect %in% non.cod]

####

# snps.dt %>%
#   filter(SNP_id %in% snpdt.obj.NoInv$SNP_id) %>%
#   filter(SNP_id %in% annotation.non.cod$SNP_id) ->
#   snps.dt.FLT
### apply filter

#### Begin analyses
  setkey(samps, locality)
  f3.o.au = foreach(test_pop = c(1:dim(tests)[1]), .combine = "rbind", .errorhandling = "remove")%do%{
       #test_pop = 4
       tests[test_pop]

      ### Reset filter
        seqResetFilter(genofile)

      ### Get data
        seqSetFilter(genofile,
                    sample.id=samps[J(t(tests[test_pop])[,1])]$sampleId,
                    variant.id=snps.dt$variant.id)

        ad <- seqGetData(genofile, "annotation/format/AD")$data
        dp <- seqGetData(genofile, "annotation/format/DP")
        #sampleids <- seqGetData(genofile, "sample.id")

        ad.dt <- as.data.table(t(ad))
        setnames(ad.dt, names(ad.dt), seqGetData(genofile, "sample.id"))
        ad.dt[,variant.id:=seqGetData(genofile, "variant.id")]

        dp.dt <- as.data.table(t(dp))
        setnames(dp.dt, names(dp.dt), seqGetData(genofile, "sample.id"))
        dp.dt[,variant.id:=seqGetData(genofile, "variant.id")]

        ad.dt.l <- melt(ad.dt, "variant.id", value.name="ad", variable.name="sampleId")
        dp.dt.l <- melt(dp.dt, "variant.id", value.name="dp", variable.name="sampleId")
        setkey(ad.dt.l, variant.id, sampleId)
        setkey(dp.dt.l, variant.id, sampleId)

        dat <- merge(ad.dt.l, dp.dt.l)

        rm(ad, dp, ad.dt, dp.dt, ad.dt.l, dp.dt.l)

        dat <- merge(dat, samps[,c("sampleId", "locality"), with=F], by="sampleId")
        dat.sum <- dat[,list(ad=sum(ad), dp=sum(dp), .N), list(variant.id, locality)]
        dat.sum.ag <- dat.sum[,list(ad=sum(ad), dp=sum(dp)), list(variant.id)]
        setkey(dat.sum, variant.id)
        dat.sum <- dat.sum[J(dat.sum.ag[ad>0 & ad!=dp & dp!=0]$variant.id)]

      ### define variables
        dat.ad.w <- as.matrix(dcast(data=dat.sum, variant.id~locality, value.var="ad"))[,-1]
        dat.dp.w <- as.matrix(dcast(data=dat.sum, variant.id~locality, value.var="dp"))[,-1]

        poolsizes <- samps[J(t(tests[test_pop])[,1])][,list(nFlies=sum(nFlies)), locality]
        setkey(poolsizes, locality)
        poolsizes <- poolsizes[J(colnames(dat.ad.w))]$nFlies

        setkey(snps.dt, variant.id)
        snp.info <- snps.dt[J(dcast(data=dat.sum, variant.id~locality, value.var="ad")$variant.id), c("chr", "pos", "ref_allele", "alt_allele"), with=F]

      ### make PoolSNP object
        pool <- new("pooldata",
                     npools=dim(dat.ad.w)[2], #### Rows = Number of pools
                     nsnp=dim(dat.ad.w)[1], ### Columns = Number of SNPs
                     refallele.readcount=dat.ad.w,
                     readcoverage=dat.dp.w,
                     poolsizes=poolsizes,
                     poolnames = colnames(dat.ad.w),
                     snp.info = snp.info)

        pool.subset<-pooldata.subset(pool, min.maf=0.05,
                                    min.cov.per.pool = 5,
                                    verbose=TRUE)

      ### f3 stats
        fstats.dat <- compute.fstats(pool.subset, computeDstat=F, nsnp.per.bjack.block = 1000)

        f3stat <- as.data.table(fstats.dat@f3star.values)
        f3stat[,focal:=tstrsplit(rownames(fstats.dat@f3star.values), ";")[[1]]]
        f3stat <- f3stat[focal==tests[test_pop]$focal]
        f3stat[,euro_parent:=tests[test_pop]$Var2]
        f3stat[,afro_parent:=tests[test_pop]$Var3]
        setnames(f3stat, names(f3stat)[1:4], c("f3", "bjack_mean", "bjack_se", "Z"))
        f3stat

    }

### save
  root = "/scratch/aob2x/f3_admix_samps/March1_2024.admix/"

  save(f3.o.au, file=paste(root, jobId, ".admix", ".Rdata", sep ="" ))
