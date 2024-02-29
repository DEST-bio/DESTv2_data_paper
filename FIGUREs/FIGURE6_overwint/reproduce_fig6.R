## Reproduce Fig 6
library(tidyverse)
library(magrittr)
library(foreach)
library(data.table)
library(gdata)
library(foreach)
library(nasapower)
library(sp)
library(lubridate)
library(stringr)
library(car)
#load geoanalysis packages
library(rnaturalearth)
library(rnaturalearthdata)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")
###
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE6_overwint/data_reproduction/ti.ob.1y.Rdata")

setDT(ti.ob.1y)
ti.ob.1y %<>%
  filter(pop1 != "Providence")

samps %>%
  group_by(pop1=city) %>%
  summarize(lat.m = mean(lat),
            long.m = mean(long)) -> lats

left_join(ti.ob.1y, lats) -> ti.ob.1y

ti.ob.1y %>%
  ggplot(aes(
    x=lat.m,
    y=logit(abs(ti.ob.1y$FST)),
    fill=T.mean,
  )) +
  geom_point(size = 3, shape = 21) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  scale_fill_gradient2(low= "steelblue", 
                       high="firebrick4", 
                       midpoint = 12) ->
  lat.fst.plot
ggsave(lat.fst.plot, 
       file = "lat.fst.plot.pdf", 
       w = 5, h = 3)
#### Add Map
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1
  ) + theme_classic() -> base_world


ti.ob.1y %>%
  group_by(pop1, lat.m, long.m) %>%
  summarise(FST_m = mean(FST)) -> FST_sum

base_world + 
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-140, 41.00), ylim = c(20.00, 68.00), expand = FALSE) + 
  geom_point(
    data = FST_sum,
    aes(x=long.m,
        y=lat.m,
        fill = FST_m), size = 3.5, shape = 21
  ) + scale_fill_gradient2(mid = "lightblue", high = "darkblue", midpoint = 0) ->
  map.overfst
ggsave(map.overfst, 
       file = "map.overfst.pdf", 
       w = 9, h = 3)

### Permute findings
perms_overw = 
  foreach(i=1:500, .combine = "rbind")%do%{
    data.frame(
      i = i,
      cor=cor.test(logit(abs(ti.ob.1y$FST)), sample(ti.ob.1y$lat.m) )$estimate 
    )
  }

rbind(
  data.frame(i = 0, cor = cor.test(logit(abs(ti.ob.1y$FST)), ti.ob.1y$lat.m )$estimate),
  perms_overw) %>% 
  mutate(perm.stat = case_when(i == 0 ~ "real",
                               TRUE ~ "perm") ) -> perm.test
## plot permutations
ggplot() +
  geom_density(data = filter(perm.test, perm.stat == "perm"),
               aes(cor), fill = "grey"
  ) + 
  geom_vline(data = filter(perm.test, perm.stat == "real"),
             aes(xintercept=cor), color = "red") +
  theme_bw() ->
  perm.plots
ggsave(perm.plots, file = "perm.plots.pdf", w = 4, h = 4)


###
## Yesiloz
####
ti.ob.1y %>%
  filter(pop1 == "Yesiloz") %>%
  mutate(year_l = paste(year1,year2)) %>%
  ggplot(aes(
    x=T.mean,
    y=logit(abs(.$FST)),
    fill=as.numeric(T.mean)) 
  ) +
  geom_point(size = 3, shape = 21) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_wrap(~year_l=="2020 2021", scale = "free_x") +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot.temp
ggsave(lat.fst.plot.temp, 
       file = "Yesiloz.lat.fst.plot.temp.pdf", 
       w = 6, h = 3)

ti.ob.1y %>%
  mutate(year_l = paste(year1,year2)) %>%
  filter(pop1 == "Yesiloz") -> Tur.1y
setDT(Tur.1y)
cor.test(Tur.1y[year_l=="2020 2021"]$T.mean, 
         Tur.1y[year_l=="2020 2021"]$FST)
cor.test(Tur.1y[year_l!="2020 2021"]$T.mean, 
         Tur.1y[year_l!="2020 2021"]$FST)

####
####
ti.ob.1y %>%
  filter(pop1 != "Yesiloz") %>%
  mutate(year_l = paste(year1,year2)) %>%
  ggplot(aes(
    x=T.mean,
    y=logit(abs(.$FST)),
    fill=as.numeric(T.mean)) 
  ) +
  geom_point(size = 3, shape = 21) +
  geom_smooth(method = "lm") +
  theme_bw() +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  Notur.lat.fst.plot.temp
ggsave(Notur.lat.fst.plot.temp, 
       file = "Notur.lat.fst.plot.temp.pdf", 
       w = 6, h = 3)

ti.ob.1y %>%
  filter(pop1 != "Yesiloz") -> NoTur.1y
cor.test(NoTur.1y$T.mean, NoTur.1y$FST)