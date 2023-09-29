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

dat.o <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))

### open metadata
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
#####
#grep( "SIM" , samps$sampleId) -> sim.pos
#grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

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
                     admix.set = "Australia"
                       ) %>%
              mutate(source_pop = rownames(.)) ->
              o
            
            return(o)
             }
     oa.g %>%
       left_join(samps) -> oa.g.m
          
  save(oa.g.m,
       file = 
         paste("linear_admix/Australia/AU.Admix", eu_pop, "Rdata", sep ="." ))     

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
                            admix.set = "N.America"
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
         paste("linear_admix/N.America/NA.Admix", eu_pop, "Rdata", sep ="." ))     
  
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
                                admix.set = "S.America"
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
         paste("linear_admix/S.America/SA.Admix", eu_pop, "Rdata", sep ="." ))     
  