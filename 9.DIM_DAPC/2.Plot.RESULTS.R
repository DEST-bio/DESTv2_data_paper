
library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(zoo)
library(adegenet)
library(reshape2)
library("geosphere")
library(gmodels)
library(ggforce)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

###
###
system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

samps <- vroom("dest_v2.samps_25Feb2023.csv")

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps

###
###
load("all.preds.Dest1.Rdata")
load("hav.dist.obj.top.pred.Rdata")
load("all.preds.Dest1.city.Rdata")
load("hav.dist.obj.top.pred.city.Rdata")
###



#### PLOT performances
rbind(mutate(hav.dist.obj.top.pred.city, model = "city|GIM 1.0"),
      mutate(hav.dist.obj.top.pred, model = "state/region|GIM 1.0")) %>%
  ggplot(aes(
    x=continent,
    y=hav_d,
    color = model
  )) +
  geom_boxplot(width = 0.8) +
  ylab(expression(italic(d)["hav"])) + 
  theme_bw() +
  coord_flip() +
  scale_y_continuous(trans='log10') ->  DAPC.v1

ggsave(DAPC.v1, file = "DAPC.v1.pdf", w= 5, h = 4)

#### PLOT
hav.dist.obj.top.pred.city %>%
  left_join(samps, by = "sampleId" ) -> pred.samps.city

samps%>%
  filter(set %in% c("dgn", "DrosEU", "DrosRTEC") ) %>%
  group_by(country,city) %>%
  summarize(a.lon = mean(long),
            a.lat = mean(lat),
            N = n()
              ) -> city.coords


CountofN_prox =
  foreach(i = 1:dim(pred.samps.city)[1], .combine = "rbind")%do%{
    pred.samps.city[i,] -> tmp
    
    city.coords %>% 
      group_by(country, city) %>%
      mutate(hav_d = distHaversine(
        matrix(c(tmp$r.long, a.lon, tmp$r.lat, a.lat),nrow = 2))/1000) %>% 
      filter(hav_d < 400) %>%
      .$N %>% sum() -> N.tmp400
    
    #city.coords %>% 
    #  group_by(country, city) %>%
    #  mutate(hav_d = distHaversine(
    #    matrix(c(tmp$r.long, a.lon, tmp$r.lat, a.lat),nrow = 2))/1000) %>% 
    #  filter(hav_d < 100) %>%
    #  .$N %>% sum() -> N.tmp100
    
    return(data.frame(tmp, N400 = N.tmp400
                      #N100 = N.tmp100
                      ))
  }

CountofN_prox %<>%
  left_join(samps, by = "sampleId" )



summary(lm(hav_d ~ N400+continent, data = CountofN_prox))
anova(lm(hav_d ~ N400+continent, data = CountofN_prox))

filter(CountofN_prox, province.y %in% 
         c("Rhode Island", "Texas", "Virginia")) %>%
  group_by(State = province.y) %>%
  summarize(N400 = mean(N400),
            hav_d =mean(hav_d),
            continent = "North_America") ->
  NAMe.examples

CountofN_prox %>%
  filter(continent %in% c("Europe", "North_America") ) %>% 
  ggplot(aes(
    y=hav_d,
    x=N400,
    #fill=variable
  )) + 
  geom_jitter(shape = 21, size = 3.1, fill = "grey") +
  geom_smooth(method = "lm", color = "black") +
  geom_point(data = NAMe.examples,
             aes(fill = State), shape = 23, size = 4.5) +
  theme_bw() +
  #theme(legend.pos = "none") +
  ylab(expression(italic(d)["hav"])) + 
  xlab("Training samples within 400 Km") +
  #scale_y_continuous(trans='log10') +
  facet_grid(continent~.) ->
  res.N400
ggsave(res.N400, file = "res.N400.pdf", w=4.5, h=4)

