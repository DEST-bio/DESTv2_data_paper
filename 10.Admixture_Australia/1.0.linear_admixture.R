# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

# module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1]) ## 1-347


### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  
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

    ### AUS
     AUS_samps = samps[continent=="Oceania"]$sampleId
     
     oa.g = foreach(au.i=AUS_samps, 
             .combine = "rbind", 
             .errorhandling = "remove")%do%{
               
               message(au.i)
          #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
            ###au.i=AUS_samps[1]
            test_pop = au.i 
            
            seqSetFilter(genofile, 
                         variant.id=sample(snps.dt[nAlleles==2]$variant.id, 1000),
                         sample.id=c(af_pop, eu_pop, test_pop))
            
            dat<- t(seqGetData(genofile, "annotation/format/FREQ")$data)
            dat <- as.data.table(dat)
            
            setnames(dat, names(dat), seqGetData(genofile, "sample.id"))
            
            setnames(dat, names(dat), c("aus", "eu", "af"))
             
            summary(lm(aus~0+eu+af, dat)) -> mod.out
            
            as.data.frame(mod.out$coefficients) %>%
              mutate(parent_eu = eu_pop,
                     parent_af = af_pop,
                     sample = test_pop
                       ) %>%
              mutate(source_pop = rownames(.)) ->
              o
            
            return(o)
             }
     
     oa.g %>%
       mutate(sampleId = sample) %>%
       left_join(samps) -> oa.g.m
          
  save(oa.g.m,
       file = 
         paste("linear_admix/AUS_Admix_EU_par", eu_pop, "Rdata", sep ="." ))     
