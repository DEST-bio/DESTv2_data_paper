#### Collect and plot admixture proportions 
library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(foreach)
require(gmodels)

meta <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")

files <- system(paste("ls moments_output_NoMigration" ), intern = T)

inf <- foreach(i=files, .combine = "rbind", .errorhandling = "remove")%do%{
  tmp <- fread(paste("moments_output_NoMigration/",i, sep = ""))
}
inf %<>% 
  mutate(model = "noMig")

files2 <- system(paste("ls moments_output_wMigration" ), intern = T)

inf2 <- foreach(i=files2, .combine = "rbind", .errorhandling = "remove")%do%{
  tmp <- fread(paste("moments_output_wMigration/",i, sep = ""))
}
inf2 %<>% 
  mutate(model = "Mig")

rbind(inf, inf2) ->
  datin_raw

datin_raw %>%
  group_by(model) %>%
  summarise(AIC = mean(AIC))

datin_raw %>%
  ggplot(aes(
    x=model,
    y=AIC
  )) + geom_boxplot() ->
  AIC_box
ggsave(AIC_box, file = "AIC_box.pdf")



datin_raw %>%
  mutate(sampleId = V1) %>%
  group_by(sampleId, model) %>%
  summarise(ancestry = ci(admix_prop)[1],
            uci = ci(admix_prop)[2],
            lci = ci(admix_prop)[3]) %>%
  left_join(meta) ->
  o

save(o, file = "../MOMENTS.AdmixtureProp.EuropeanProb.Rdata")

o %>%
  ggplot(aes(
    x=lat,
    y=1-ancestry,
    color=model
  )) + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~continent, scales = "free_x") +
  theme_bw() -> momemts_ests

ggsave(momemts_ests, file = "momemts_ests.pdf",
       w= 9, h = 3)

#####
lndat <- "/gpfs2/scratch/jcnunez/DEST2.0_analysis/f3_revist/linear.admix.dat.Rdata"

lindat <- get(load(lndat))
lindat %>% filter(filter=="noINV") %>%
  filter(source_pop == "AFRICA") %>%
  ungroup() %>%
  dplyr::select(sampleId, linearEST = mean.est)->
  lndat.f

o %>% left_join(lndat.f) ->
  o.j

o.j %>%
  ggplot(aes(
    x=1-ancestry,
    y=linearEST
  )) + geom_point() +
  geom_abline(slope = 1) +
  facet_grid(~continent, scales = "free_x") +
  theme_bw() -> corr_methods

ggsave(corr_methods, file = "corr_methods.pdf",
       w= 9, h = 3)

##### Final dataste
save(o.j, file = "ancestry_Estimates_linear_moments.Rdata")



  
  
  
