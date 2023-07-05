# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  
  
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggExtra)
  library(foreach)
  
### open GDS
aus <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/f3_Admix/*", intern = T)
nam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/f3_Admix/*", intern = T)
sam <- system("ls /scratch/yey2sn/DEST2_analysis/admix_samps/f3_Admix/*", intern = T)

samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))


aus.samps = 
  foreach(fil = aus, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } 

nam.samps = 
  foreach(fil = nam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } 

sam.samps = 
  foreach(fil = sam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } 


rbind(
aus.samps,
nam.samps,
sam.samps) ->
  f3.admix.dat

save(f3.admix.dat, file = "f3.admix.dat.Rdata")
####
f3.admix.dat %>%
  filter(admix.set == "Australia") %>%
  lm(f3 ~ lat, data = .) %>% summary()
f3.admix.dat %>%
  filter(admix.set == "S.America") %>%
  lm(f3 ~ lat, data = .) %>% summary()
f3.admix.dat %>%
  filter(admix.set == "N.America") %>%
  lm(f3 ~ lat, data = .) %>% summary()

####
f3.admix.dat %>% 
  left_join(. , dplyr::select(samps,parent_eu=sampleId, parent_lat=lat,
                              parent_long=long, parent_country = country)) ->
  f3.admix.dat.p

####
f3.admix.dat.p %>%
  group_by(parent_long, parent_lat, admix.set) %>%
  summarize(m.f3 = sum(f3 < 0),
            N = n()) %>% 
  filter(m.f3/N != 0) ->
  f3.admix.dat.p.ag
  
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

  ggplot(data = world) +
  geom_sf(fill= "grey70") +
    coord_sf(xlim = c(-12.50, 41.00), ylim = c(33.00, 69.00), expand = FALSE) + 
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "white")) +
geom_jitter(data = f3.admix.dat.p.ag , 
           aes(
  x=parent_long,
  y=parent_lat,
  #size = m.f3/N,
  fill = m.f3/N
), alpha = 0.9, size = 2.1, shape = 21, color = "white") + 
  scale_fill_gradient2(low = "lightblue", mid = "blue", high = "darkblue", midpoint = 0.45) +
  facet_grid(~admix.set)->
  plot.f3admix

ggsave(plot.f3admix, file = "plot.f3admix.png", w = 9, h  =3)



