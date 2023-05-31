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

### metadata
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"
samps <- fread(meta_git)
setDT(samps)

### open GDS
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

####
snpdt.obj <- get(load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt.Rdata"))
setDT(snpdt.obj)
snpdt.obj %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
snpdt.obj.NoInv = snpdt.obj[invName == "none"]

annotation <- get(load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_Cville_SNP_Metadata.Rdata"))
setDT(annotation)
names(annotation)[1:2] = c("chr","pos")
non.cod = c("intergenic_region","intron_variant","upstream_gene_variant","downstream_gene_variant")
annotation.non.cod = annotation[effect %in% non.cod]
####

snps.dt %>% 
  filter(SNP_id %in% snpdt.obj.NoInv$SNP_id) %>%
  filter(SNP_id %in% annotation.non.cod$SNP_id) ->
  snps.dt.FLT
### apply filter

#### Begin analyses
#### 
#### 

  af_pop <- samps[continent=="Africa"][nFlies>100]$sampleId
  eu_pop <- samps[continent=="Europe"]$sampleId[k]
  
  children_samps = samps[continent=="North_America"]$sampleId

#### ----> Run loop
f3.o.au = foreach(ch.i = children_samps, 
               .combine = "rbind", 
               .errorhandling = "remove")%do%{
                 
                 #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                 ###au.i=AUS_samps[1]
                 ch.i = children_samps[1]
                 test_pop = ch.i
                 message(ch.i)
                 
                 ### Reset filter
                 seqResetFilter(genofile)
                 
                 ### Get data
                 seqSetFilter(genofile, sample.id=c(af_pop,eu_pop,test_pop),
                              variant.id=snps.dt.FLT$variant.id)
                 
                 ####
                 ad <- seqGetData(genofile, "annotation/format/AD")$data
                 dp <- seqGetData(genofile, "annotation/format/DP")
                 sampleids <- seqGetData(genofile, "sample.id")
                 ####
                 ####
                 colnames(ad) <- paste(seqGetData(genofile, "chromosome"), 
                                        seqGetData(genofile, "position") , 
                                        sep="_")
                 colnames(dp) <- paste(seqGetData(genofile, "chromosome"), 
                                       seqGetData(genofile, "position") , 
                                       sep="_")
                 ####
                 ####
                 #### ----> Some cleaning up :)
                 apply(ad, 2,  function(x) { sum( x == 0)}  ) %>% 
                   data.frame(fixed = .) %>% 
                   filter(fixed == 1) %>% 
                   row.names(.) -> poly.sites
                 
                 ####
                 ####
                 ad.flt = ad[,poly.sites]
                 dp.flt = dp[,poly.sites]
                
                 snps.dt.FLT %>%
                   filter(SNP_id %in% poly.sites) ->
                   snps.dt.ply
                 ####
                 ####
                
                 poolsizes = c(
                   filter(samps, sampleId %in% sampleids)$nFlies
                 )
                 
                 ####
                 ####
                 ####
                 ####
                 pool <- new("pooldata",
                             npools=dim(ad.flt)[1], #### Rows = Number of pools
                             nsnp=dim(dp.flt)[2], ### Columns = Number of SNPs
                             refallele.readcount=t(ad.flt),
                             readcoverage=t(dp.flt),
                             poolsizes=poolsizes,
                             poolnames = sampleids,
                             snp.info = snps.dt.ply[,c("chr","pos","ref_allele","alt_allele")])
                 
                 pool.subset<-pooldata.subset(pool, min.maf=0.05,
                                              min.cov.per.pool = 10,
                                              verbose=TRUE)
                 
                 fstats.dat <- 
                   compute.fstats(
                     pool.subset,
                     computeDstat=TRUE,
                     nsnp.per.bjack.block = 1000
                   )
                 
                 fstats.dat@f3star.values %>%
                   as.data.frame()  ->
                   f3stat
                 f3stat
                 
                 target.pop.col = grep(paste(test_pop,";", sep = ""), rownames(f3stat))
                 data.frame(f3 = f3stat[target.pop.col,]) %>%
                   mutate(set.tree = rownames(f3stat)[target.pop.col] ) %>%
                   mutate(parent_eu = eu_pop,
                          parent_af = af_pop,
                          sampleId = test_pop)->
                   o
                 
                 return(o)
               }

####
NAM_samps = samps[continent=="North_America"]$sampleId

f3.o.Nam = foreach(NAM.i=NAM_samps, 
                   .combine = "rbind", 
                   .errorhandling = "remove")%do%{
                     
                     message(NAM.i)
                     #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                     ###au.i=AUS_samps[1]
                     test_pop = NAM.i 
                    
                    poolsizes = c(
                      filter(samps, sampleId == rownames(dp)[1])$nFlies,
                      filter(samps, sampleId == rownames(dp)[2])$nFlies,
                      filter(samps, sampleId == rownames(dp)[3])$nFlies
                    )
                    
                    
                    ad.tmp = ad[ c(af_pop, eu_pop, test_pop), sampled.SNPs]
                    
                    dp.tmp = dp[ c(af_pop, eu_pop, test_pop), sampled.SNPs]
                    
                    pool <- new("pooldata",
                                npools=dim(ad.tmp)[1], #### Rows = Number of pools
                                nsnp=dim(ad.tmp)[2], ### Columns = Number of SNPs
                                refallele.readcount=t(ad.tmp),
                                readcoverage=t(dp.tmp),
                                poolsizes=poolsizes * 2,
                                poolnames = rownames(ad.tmp),
                                snp.info = snp.info)
                    
                    fstats.dat <- 
                      compute.fstats(
                        pool,
                        #nsnp.per.bjack.block = 100,
                        computeDstat = TRUE,
                        verbose = TRUE
                      )
                    
                    fstats.dat@f3star.values %>%
                      as.data.frame()  ->
                      f3stat
                    
                    target.pop.col = grep(paste(test_pop,";", sep = ""), rownames(f3stat))
                    data.frame(f3 = f3stat[target.pop.col,]) %>%
                      mutate(set.tree = rownames(f3stat)[target.pop.col] ) %>%
                      mutate(parent_eu = eu_pop,
                             parent_af = af_pop,
                             sampleId = test_pop)->
                      o
                    
                    return(o)
                  }

#####
SAM_samps = samps[continent=="South_America"]$sampleId

f3.o.Sam = foreach(SAM.i=SAM_samps, 
                   .combine = "rbind", 
                   .errorhandling = "remove")%do%{
                     
                     message(SAM.i)
                     #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                     ###au.i=AUS_samps[1]
                     test_pop = SAM.i 
                     
                     poolsizes = c(
                       filter(samps, sampleId == rownames(dp)[1])$nFlies,
                       filter(samps, sampleId == rownames(dp)[2])$nFlies,
                       filter(samps, sampleId == rownames(dp)[3])$nFlies
                     )
                     
                     
                     ad.tmp = ad[ c(af_pop, eu_pop, test_pop), sampled.SNPs]
                     
                     dp.tmp = dp[ c(af_pop, eu_pop, test_pop), sampled.SNPs]
                     
                     pool <- new("pooldata",
                                 npools=dim(ad.tmp)[1], #### Rows = Number of pools
                                 nsnp=dim(ad.tmp)[2], ### Columns = Number of SNPs
                                 refallele.readcount=t(ad.tmp),
                                 readcoverage=t(dp.tmp),
                                 poolsizes=poolsizes * 2,
                                 poolnames = rownames(ad.tmp),
                                 snp.info = snp.info)
                     
                     fstats.dat <- 
                       compute.fstats(
                         pool,
                         #nsnp.per.bjack.block = 100,
                         computeDstat = TRUE,
                         verbose = TRUE
                       )
                     
                     fstats.dat@f3star.values %>%
                       as.data.frame()  ->
                       f3stat
                     
                     target.pop.col = grep(paste(test_pop,";", sep = ""), rownames(f3stat))
                     data.frame(f3 = f3stat[target.pop.col,]) %>%
                       mutate(set.tree = rownames(f3stat)[target.pop.col] ) %>%
                       mutate(parent_eu = eu_pop,
                              parent_af = af_pop,
                              sampleId = test_pop)->
                       o
                     
                     return(o)
                   }

####

rbind(
mutate(f3.o.au, admix.set = "Australia"),
mutate(f3.o.Nam, admix.set = "N.America"),
mutate(f3.o.Sam, admix.set = "S.America")) %>%
  left_join(samps) -> out.tmp.merged


save(out.tmp.merged,
     file = 
       paste("/scratch/yey2sn/DEST2_analysis/admix_samps/f3_Admix/", eu_pop, ".Rdata", sep ="" ))
