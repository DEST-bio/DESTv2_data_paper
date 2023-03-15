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

dp <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/DPmatrix.flt.Rdata"))
ad <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/ADmatrix.flt.Rdata"))
colnames(ad)[sample(dim(ad)[2],10000)] -> sampled.SNPs
snp.info = sampled.SNPs %>%
  data.frame(SNP_id=.) %>%
  separate(SNP_id, into = c("chr","pos","id")) %>%
  dplyr::select(chr, pos)
  
### open metadata
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
#####
#grep( "SIM" , samps$sampleId) -> sim.pos
#grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

### define test pops
af_pop <- samps[continent=="Africa"][nFlies>100]$sampleId
eu_pop <- samps[continent=="Europe"]$sampleId[k]

AUS_samps = samps[continent=="Oceania"]$sampleId

f3.o.au = foreach(au.i=AUS_samps, 
               .combine = "rbind", 
               .errorhandling = "remove")%do%{
                 
                 message(au.i)
                 #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                 ###au.i=AUS_samps[1]
                 test_pop = au.i 
                 
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
       paste("/yey2sn/DEST2_analysis/admix_samps/f3_Admix/", eu_pop, ".Rdata", sep ="" ))     

