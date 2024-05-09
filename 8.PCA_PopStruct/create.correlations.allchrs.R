### Extra PCA analyses
library(tidyverse)
library(foreach)
library(reshape2)

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE3_PCA/data_for_reproduction/PCA.results.df.Rdata")

PCA.results.df$case %>% table

cors.pcs.chrs=
foreach( ch = c("2L", "2R", "3L", "3R"), 
         .combine = "rbind")%do%{
           
           message(ch)
           PCA.results.df %>%
             filter(case == ch) ->
             tmp
           
           foreach(cont = c("Europe", "North_America"), 
                   .combine = "rbind")%do%{
             message(cont)
             
             tmp %>%
               filter(continent == cont) ->
               tmp2
             
            la1= cor.test(tmp2$Dim.1, tmp2$lat)
            la2= cor.test(tmp2$Dim.2, tmp2$lat)
            la3= cor.test(tmp2$Dim.3, tmp2$lat)
             
            lo1= cor.test(tmp2$Dim.1, tmp2$long)
            lo2= cor.test(tmp2$Dim.2, tmp2$long)
            lo3= cor.test(tmp2$Dim.3, tmp2$long)
           
            data.frame(
              ch,
              cont,
              la1_e = la1$estimate,
              la1_lci = la1$conf.int[1],
              la1_uci = la1$conf.int[2],
              
              la2_e = la2$estimate,
              la2_lci = la2$conf.int[1],
              la2_uci = la2$conf.int[2],
              
              la3_e = la3$estimate,
              la3_lci = la3$conf.int[1],
              la3_uci = la3$conf.int[2],
              
              lo1_e = lo1$estimate,
              lo1_lci = lo1$conf.int[1],
              lo1_uci = lo1$conf.int[2],
              
              lo2_e = lo2$estimate,
              lo2_lci = lo2$conf.int[1],
              lo2_uci = lo2$conf.int[2],
              
              lo3_e = lo3$estimate,
              lo3_lci = lo3$conf.int[1],
              lo3_uci = lo3$conf.int[2]
              )
            
           }}

cors.pcs.chrs %>%
  melt(id = c("ch", "cont")) %>%
  separate(variable, into = c("var", "est"), sep = "_") %>% 
  dcast(ch+cont+var~est, value.var = "value") %>%
  separate(var, into = c("vat","pc"), 2) %>%
  ggplot(aes(
    x=paste(vat,pc),
    color=pc,
    y=e,
    ymin=lci,
    ymax=uci,
  )) +
  geom_hline(yintercept =  0) +
  geom_errorbar(width = 0.1,position=position_dodge(width=0.5) ) + 
  geom_point(position=position_dodge(width=0.5)) + 
  facet_grid(cont~ch) 

