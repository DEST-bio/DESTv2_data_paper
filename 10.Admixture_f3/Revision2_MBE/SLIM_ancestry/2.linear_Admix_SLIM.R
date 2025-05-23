library(foreach)
library(data.table)
library(gmodels)
library(foreach)
library(tidyverse)
library(magrittr)


args = commandArgs(trailingOnly=TRUE)
job=as.numeric(args[1])

files <- system(paste("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/results/AFR.",job,".*", sep = ""),
                intern = T)

afr_mks =
  foreach(i = files,
          .combine = "rbind")%do%{
            
            replicate=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][2]
            pop=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][4]
            set=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][1]
            
            tmp <- fread(i) %>% t()
            data.frame(AF=tmp) %>%
              mutate(pop=pop) %>%
              mutate(set=set) %>%
              mutate(replicate=replicate)
          }


afr_mks %>%
  group_by(sampleId=pop, replicate) %>%
  summarise(ancestry = ci(AF)[1],
            lci = ci(AF)[2],
            uci = ci(AF)[3]
            ) ->
  afr_mks.ag

afr_mks.ag$sampleId = as.numeric(afr_mks.ag$sampleId)

files2 <- system(paste("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/results/Agnostic.",job,".*", sep = ""),
                intern = T)
o =
foreach(i = c(files2),
        .combine = "rbind")%do%{
    
          replicate=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][2]
          pop=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][4]
          set=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][1]
          
          tmp <- fread(i) %>% t()
          data.frame(AF=tmp) %>%
            mutate(pop=pop) %>%
            mutate(set=set) %>%
            mutate(replicate=replicate)
      
        }

o %>%
  filter(pop == 0) %>%
  mutate(pop = "AFRICA")->
  AFR

o %>%
  filter(pop == 1)%>%
  mutate(pop = "EUROPE")->
  EUR



SIM.oa.g = foreach(anchor=2:(length(files)-1), 
                   .combine = "rbind")%do%{
                     
                     message(anchor)
                     o %>%
                       filter(pop == anchor,
                              set == "Agnostic") %>%
                       mutate(pop="ANCHOR")->
                       ANCHOR
                     
                     data.frame(
                     anchor=ANCHOR$AF,
                     europe=EUR$AF,
                     africa=AFR$AF) ->
                       merged_AF
                     
                     mAF = apply(merged_AF, 1, function(x) mean(x)  )
                      
                     merged_AF %<>%
                       mutate(mAF = mAF)
                     
                     
                     merged_AF %>%
                       filter(mAF > 0.05) ->
                       dat
                    
                     summary(lm(anchor~0+europe+africa, data = dat)) -> mod.out
                     
                     as.data.frame(mod.out$coefficients) %>%
                       mutate(parent_1 = "0",
                              parent_2 = "1",
                              sampleId = anchor,
                              admix.set = "Simulations"
                       ) %>%
                       mutate(source_pop = rownames(.)) ->
                       o2
                     
                     return(o2)
                   }


SIM.oa.g %>%
  filter(source_pop == "africa") %>%
  left_join( afr_mks.ag) %>%
  mutate(rep=job)->
  coeffs

system("mkdir ag.results")
save(coeffs,
     file = paste("ag.results/",job,".ag.ancestry.Rdata", sep = "" ))

#coeffs %>%
#  ggplot(aes(
#    x=sampleId,
#    
#  )) + geom_point(aes(y=Estimate)) +
#  geom_smooth(aes(y=Estimate), method = "lm") +
#  geom_point(aes(y=ancestry), color = "red") +
#  geom_errorbar(aes(ymin=lci, ymax=uci), color = "red", width = 0.1) +
#  ggtitle("SLiM simulation", subtitle = "Red dots are the 'true' ancestry")
#


