# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

# module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1]) ## 1-347
filto=as.character(args[2]) # chr2L, chr2R, chr3L, chr3R, noINV, 2Lt, 2RNS, 3RK, 3RK_3RP_3RMo, silent


### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  
### open GDS

dat.o <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))

### open metadata
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
#####
#grep( "SIM" , samps$sampleId) -> sim.pos
#grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos
####

snpdt.obj <- get(load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt.Rdata"))
setDT(snpdt.obj)
snpdt.obj %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))

annotation <- get(load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_Cville_SNP_Metadata.Rdata"))
setDT(annotation)

annotation %>% filter(!Putative_impact %in% c("LOW","MEDIUM","HIGH")) %>% .$SNP_id -> func.ids

### Annotation Filters

# chr2L, chr2R, chr3L, chr3R, noINV, 2Lt, 2RNS, 3RK, 3RK_3RP_3RMo, silent

if(filto == "chr2L"){
  filt.i = snpdt.obj[chr == "2L"]
} else if(filto == "chr2R"){
  filt.i = snpdt.obj[chr == "2R"]
} else if(filto == "chr3L"){
  filt.i = snpdt.obj[chr == "3L"]
} else if(filto == "chr3R"){
  filt.i = snpdt.obj[chr == "3R"]
} else if(filto == "noINV"){
  filt.i = snpdt.obj[invName == "none"]
} else if(filto == "2Lt"){
  filt.i = snpdt.obj[invName == "2Lt"]
} else if(filto == "2RNS"){
  filt.i = snpdt.obj[invName == "2RNS"]
} else if(filto == "3RK"){
  filt.i = snpdt.obj[invName == "3RK"]
} else if(filto == "3RK_3RP_3RMo"){
  filt.i = snpdt.obj[invName == "3RK;3RP;3RMo"]
} else if(filto == "silent"){
  snpdt.obj %>%
    filter(!SNP_id %in% func.ids) ->
    filt.i
} else if(filto == "all"){
  message("using all")
}

####
if(filto == "all"){
  dat.o = dat.o
} else if(filto != "all"){
  
  data.frame(snpnam = colnames(dat.o))%>%
    separate(snpnam, remove = F, into = c("chr","pos","id")) ->
    tmp.dt
  tmp.dt %<>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))
  
  filter.app = tmp.dt %>% filter(SNP_id %in% filt.i$SNP_id) %>% .$snpnam
  
  dat.o = dat.o[,filter.app] 
  message(dim(dat.o)[2])
}


###
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
            
            dat.o[ c(af_pop, eu_pop, test_pop), sample(dim(dat.o)[2],5000) ] -> 
              dat.tmp
            
            dat.tmp %>%
              t %>%
              as.data.frame() ->
              dat.t
            
            setnames(dat.t, names(dat.t), c("AFRICA", "EUROPE", "TEST"))
             
            summary(lm(TEST~0+EUROPE+AFRICA, dat.t)) -> mod.out
            
            as.data.frame(mod.out$coefficients) %>%
              mutate(parent_eu = eu_pop,
                     parent_af = af_pop,
                     sampleId = test_pop,
                     admix.set = "Australia",
                     filter = filto
                       ) %>%
              mutate(source_pop = rownames(.)) ->
              o
            
            return(o)
             }
     oa.g %>%
       left_join(samps) -> oa.g.m
          
  save(oa.g.m,
       file = 
         paste("linear_admix/Australia/AU.Admix", filto, eu_pop, "Rdata", sep ="." ))     

  oa.g.m %>%
    ggplot(aes(
      x= lat,
      y = Estimate,
      color = source_pop
    )) +
      geom_smooth(method = "lm") +
    geom_point() ->
    aus.plot
  ggsave(aus.plot, file = "aus.plot.pdf")
    
  
  ### NAM.
  NAM_samps = samps[continent=="North_America"]$sampleId
  
  NAM.oa.g = foreach(NAM.i=NAM_samps, 
                 .combine = "rbind", 
                 .errorhandling = "remove")%do%{
                   
                   message(NAM.i)
                   #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                   ###au.i=AUS_samps[1]
                   test_pop = NAM.i 
                   
                   dat.o[ c(af_pop, eu_pop, test_pop), sample(dim(dat.o)[2],5000) ] -> 
                     dat.tmp
                   
                   dat.tmp %>%
                     t %>%
                     as.data.frame() ->
                     dat.t
                   
                   setnames(dat.t, names(dat.t), c("AFRICA", "EUROPE", "TEST"))
                   
                   summary(lm(TEST~0+EUROPE+AFRICA, dat.t)) -> mod.out
                   
                   as.data.frame(mod.out$coefficients) %>%
                     mutate(parent_eu = eu_pop,
                            parent_af = af_pop,
                            sampleId = test_pop,
                            admix.set = "N.America",
                            filter = filto
                     ) %>%
                     mutate(source_pop = rownames(.)) ->
                     o
                   
                   return(o)
                 }
  
  NAM.oa.g %>%
    left_join(samps) -> NAM.oa.g.m
  
  #NAM.oa.g.m %>%
  #  ggplot(aes(
  #    x= lat,
  #    y = Estimate,
  #    color = source_pop
  #  )) +
  #    geom_smooth(method = "lm") +
  #  geom_point() ->
  #  NAM.oa.g.plot
  #ggsave(NAM.oa.g.plot, file = "NAM.oa.g.plot.pdf")
  
  save(NAM.oa.g.m,
       file = 
         paste("linear_admix/N.America/NA.Admix", filto ,eu_pop, "Rdata", sep ="." ))     
  
 ## 
 ## SAM.
  SAM_samps = samps[continent=="South_America"]$sampleId
  
  SAM.oa.g = foreach(SAM.i=SAM_samps, 
                     .combine = "rbind", 
                     .errorhandling = "remove")%do%{
                       
                       message(SAM.i)
                       #test_pop <- sample(samps[continent=="Oceania"]$sampleId, 1)
                       ###au.i=AUS_samps[1]
                       test_pop = SAM.i 
                       
                       dat.o[ c(af_pop, eu_pop, test_pop), sample(dim(dat.o)[2],5000) ] -> 
                         dat.tmp
                       
                       dat.tmp %>%
                         t %>%
                         as.data.frame() ->
                         dat.t
                       
                       setnames(dat.t, names(dat.t), c("AFRICA", "EUROPE", "TEST"))
                       
                       summary(lm(TEST~0+EUROPE+AFRICA, dat.t)) -> mod.out
                       
                       as.data.frame(mod.out$coefficients) %>%
                         mutate(parent_eu = eu_pop,
                                parent_af = af_pop,
                                sampleId = test_pop,
                                admix.set = "S.America",
                                filter = filto
                         ) %>%
                         mutate(source_pop = rownames(.)) ->
                         o
                       
                       return(o)
                     }
  
  
  SAM.oa.g %>%
    left_join(samps) -> SAM.oa.g.m
  
  #SAM.oa.g.m %>%
  #  ggplot(aes(
  #    x= lat,
  #    y = Estimate,
  #    color = source_pop
  #  )) +
  #  geom_smooth(method = "lm") +
  #  geom_point() ->
  #  SAM.oa.g.plot
  #ggsave(SAM.oa.g.plot, file = "SAM.oa.g.plot.pdf")
  
  save(SAM.oa.g.m,
       file = 
         paste("linear_admix/S.America/SA.Admix", filto,  eu_pop, "Rdata", sep ="." ))     
  