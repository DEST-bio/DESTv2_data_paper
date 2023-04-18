# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  
### open GDS
aus <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/Australia/*", intern = T)
nam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/N.America/*", intern = T)
sam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/S.America/*", intern = T)

aus.samps = 
  foreach(fil = aus, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

nam.samps = 
  foreach(fil = nam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

sam.samps = 
  foreach(fil = sam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))


rbind(
aus.samps,
nam.samps,
sam.samps) ->
  linear.admix.dat

save(linear.admix.dat, file = "linear.admix.dat.Rdata")
####
linear.admix.dat %>%
  filter(admix.set == "Australia" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat, data = .) %>% summary()
linear.admix.dat %>%
  filter(admix.set == "S.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat, data = .) %>% summary()
linear.admix.dat %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat, data = .) %>% summary()



linear.admix.dat %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  arrange(lat)
####
linear.admix.dat %>%
ggplot(aes(
  x=lat,
  y=mean.est,
  ymin = mean.est - sd.est,
  ymax = mean.est + sd.est,
  color = source_pop
)) + 
  geom_smooth(method = "lm") +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2.1, shape = 21, 
             color = "black", aes(fill = source_pop)) +
  theme_bw() +
  facet_grid(~admix.set, scales = "free_x")->
  plot.admix

ggsave(plot.admix, file = "plot.au.pdf", w = 9, h  =3)



