#### DEST Analysis FST across time and across space
#### Script 1 -- find the samples used for the FST analysis.

rm(list = ls())

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(foreach)
library(gtools)
library(doParallel)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(doMC)
registerDoMC(4)
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

ingds = "/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"

##### step 1 
### filter samps by pass
#### Filter by cluster
#### Cluster 2 == EUW
#### Cluster 8 == EUE

samps %>%
filter(Recommendation == "Pass")  ->
samps.pass

#### Define comparison sets
samps.pass$locality -> locs.all

y1comps = 
foreach(i = locs.all,
        .combine = "rbind"
        )%do%{

samps.pass %>%
filter(locality == i) %>% 
            group_by(sampleId) %>%
    mutate(collectionDate = 
             as.Date(str_split(sampleId, "_")[[1]][5], 
    format = "%Y-%m-%d")) -> tmp

if(length(unique(tmp$year)) > 1 ){
  message("pass OW filter!")
  LIN = dim(tmp)[1]
  
  comp_vector.in = combinations(
    LIN,
    2, 
    v=1:LIN,
    set=TRUE, 
    repeats.allowed=FALSE) %>%
    as.data.frame() %>%
    mutate(day_diff = NA,
           pop1 = NA,
           pop2 = NA,
           samp1 = NA,
           samp2 = NA
           )
  
    for(k in 1:dim(comp_vector.in)[1]){
    
    date1=as.Date(tmp$collectionDate[comp_vector.in[k,1]])
    date2=as.Date(tmp$collectionDate[comp_vector.in[k,2]])
    
    comp_vector.in$day_diff[k] = abs(as.numeric(date1-date2))
    
    comp_vector.in$pop1[k] = tmp$city[comp_vector.in[k,1]]
    comp_vector.in$pop2[k] = tmp$city[comp_vector.in[k,2]]
    
    comp_vector.in$samp1[k] = tmp$sampleId[comp_vector.in[k,1]]
    comp_vector.in$samp2[k] = tmp$sampleId[comp_vector.in[k,2]]
    
  }
  
  comp_vector.in %>%
    group_by(samp1,samp2) %>%
    mutate(year1 = 
             year(as.Date(str_split(samp1, "_")[[1]][5], 
                          format = "%Y-%m-%d")),
           year2 = 
             year(as.Date(str_split(samp2, "_")[[1]][5], 
                          format = "%Y-%m-%d")),
    ) %>%
    mutate(ydiff = abs(year1-year2)) ->
    comp_vector.in.y
  
  out =
    comp_vector.in.y %>% filter(ydiff == 1)
  
} else{
  message("not enough years")
}

}

save(y1comps,
     file = "y1comps.guide.Rdata")

#### end

