# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  library(reshape2)
  
### open GDS
aus <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/Australia/*", intern = T)
nam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/N.America/*", intern = T)
sam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/linear_admix/S.America/*", intern = T)

aus.samps = 
  foreach(fil = aus, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

nam.samps = 
  foreach(fil = nam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

sam.samps = 
  foreach(fil = sam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))


rbind(
aus.samps,
nam.samps,
sam.samps) ->
  linear.admix.dat.filters

save(linear.admix.dat.filters, file = "linear.admix.dat.Rdata")
####
regression.coeffs = 
foreach(filter.i = unique(linear.admix.dat.filters$filter),
        .combine = "rbind")%do%{
        
linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "Australia" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary() -> o1
          
linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "S.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary()  -> o2

linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary()  -> o3

rbind(
data.frame(
mod = "Australia",
lat = o1$coefficients[2,1],
p = o1$coefficients[2,4],
filter =  filter.i),
data.frame(
mod = "S.America",
lat = o2$coefficients[2,1],
p = o2$coefficients[2,4],
filter =  filter.i),
data.frame(
mod = "N.America",
lat = o3$coefficients[2,1],
p = o3$coefficients[2,4],
filter =  filter.i))

        }

regression.coeffs %>%
  reshape2::dcast(filter~mod, value.var = "p")

linear.admix.dat %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  arrange(lat)
####
linear.admix.dat.filters %>%
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
  facet_grid(filter~admix.set, scales = "free_x")->
  plot.admix.flt

ggsave(plot.admix.flt, file = "plot.admix.flt.pdf", w = 9, h  =6)



