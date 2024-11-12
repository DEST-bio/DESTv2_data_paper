## Reproduce fig 5

library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(gmodels)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

setwd("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE7_GIMS/dat_reprod")

all.preds.Dest1 = get(load("all.preds.Dest1.Rdata"))
hav.dist.obj.top.pred = get(load("hav.dist.obj.top.pred.Rdata"))
all.preds.Dest1.city = get(load("all.preds.Dest1.city.Rdata"))
hav.dist.obj.top.pred.city = get(load("hav.dist.obj.top.pred.city.Rdata"))

### comapre
left_join(hav.dist.obj.top.pred %>%
            dplyr::select(sampleId, Region_pos=posterior, Region_havd=hav_d),
          hav.dist.obj.top.pred.city %>%
            dplyr::select(sampleId, City_pos=posterior, City_havd=hav_d)) -> DEST1.comps

DEST1.comps %>%
  ggplot(aes(
    x=log10(Region_havd),
    y=log10(City_havd)
  )) + 
  geom_point(size = 2.5, shape = 21, fill = "grey") +
  theme_bw() +
  geom_smooth(method = "lm") ->
  pos.cor
ggsave(pos.cor, file = "DEST1.havd.cor.pdf", w= 4, h = 4)


###
load("DEST2.0.model.predictions.Rdata")

### Plot Map
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1
  ) + theme_classic() -> base_world

gims.out.v2 %>%
  group_by(real.lat, real.long) %>%
  summarise(dist.to.real = median(hav_d)) -> hav_d_mean

setDT(hav_d_mean)
base_world + 
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-135, 158.00), ylim = c(-48.00, 68.00), expand = FALSE) + 
  geom_jitter(
    data = hav_d_mean,
    aes(y=real.lat,
        x=real.long,
        fill = log10(dist.to.real+1)), size = 2.9, shape = 21, alpha = 0.8
  ) + scale_fill_gradient2(midpoint = 2.1, low = "yellow", high = "red") ->
  prediction.quality
ggsave(prediction.quality, 
       file = "prediction.quality.pdf", 
       w = 12, h = 4)
####


rbind(
  mutate(hav.dist.obj.top.pred[,c("continent","hav_d")], model = "DEST 1.0"), 
  mutate(
    dplyr::select(gims.out.v2, continent=cont,hav_d )[,-1], model = "DEST 2.0")) %>%
  group_by(continent, model) %>%
  summarize(#mean.d = ci(hav_d)[1],
            #lci= ci(hav_d)[2],
            #uci= ci(hav_d)[3],
            mean.d = quantile(hav_d, 0.5),
            lci= quantile(hav_d, 0.1),
            uci= quantile(hav_d, 0.9),
            m.sd = sd(hav_d)
  ) -> mod.summaries

d1.med = median(filter(mod.summaries, model == "DEST 1.0")$mean.d)
d2.med = median(filter(mod.summaries, model == "DEST 2.0")$mean.d)

mod.summaries$lci[mod.summaries$lci < 0]=1

mod.summaries %>%
  ggplot(aes(
    x=continent,
    y=mean.d,
    ymin = lci,
    ymax = uci,
    color = model
  )) +
  geom_hline(yintercept = d1.med, 
             linetype = "dashed", color = "firebrick") +
  geom_hline(yintercept = d2.med, 
             linetype = "dashed", color = "steelblue") +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  ylab(expression(italic(d)["hav"])) + 
  theme_bw() +
  scale_color_manual(values = c("firebrick","steelblue")) +
  coord_flip() + scale_y_continuous(trans = "log10") ->  DAPC.v1.v2
#scale_y_continuous(trans='log10') 

ggsave(DAPC.v1.v2, file = "DAPC.v1.v2.model.performance.pdf", w= 3.5, h = 4)

# Assessment of distance and predictive power

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

o_50_1500 =
  foreach(DIST = seq(from =50, to = 1500, by = 100),
          .combine = "rbind")%do%{
            
            CountofN_prox =
              foreach(i = 1:dim(pred.samps.city)[1], .combine = "rbind")%do%{
                pred.samps.city[i,] -> tmp
                
                message(DIST)
                
                city.coords %>% 
                  group_by(country, city) %>%
                  mutate(hav_d = distHaversine(
                    matrix(c(tmp$r.long, a.lon, tmp$r.lat, a.lat),nrow = 2))/1000) %>% 
                  filter(hav_d < DIST) %>%
                  .$N %>% sum() -> N.tmp400
                
                return(data.frame(tmp, N400 = N.tmp400
                                  #N100 = N.tmp100
                ))
              }
            
            CountofN_prox %<>%
              left_join(samps, by = "sampleId" )
            
            
            summary(lm(hav_d ~ N400+continent, data = CountofN_prox)) -> lm.o
            
            o =
              data.frame(
                beta = lm.o$coefficients[2],
                dist = DIST)
            
            return(o)
          }

####
o_50_1500 %>%
  ggplot(
    aes(
      x=dist,
      y=beta*-1)) +
  theme_bw() +
  geom_smooth() +
  geom_point(size = 2.5, shape = 21, fill = "grey") ->
  improvement_model

ggsave(improvement_model, file = "improvement_model.pdf", w= 4, h = 4)

#### a special case
CountofN_0 =
  foreach(i = 1:dim(pred.samps.city)[1], .combine = "rbind")%do%{
    pred.samps.city[i,] -> tmp
    
    city.coords %>% 
      group_by(country, city) %>%
      mutate(hav_d = distHaversine(
        matrix(c(tmp$r.long, a.lon, tmp$r.lat, a.lat),nrow = 2))/1000) %>% 
      filter(hav_d <1) %>%
      .$N %>% sum() -> N.tmp
    
    return(data.frame(tmp, N0 = N.tmp, D = 1
                      #N100 = N.tmp100
    ))
  }

CountofN_0 %>%
  filter(N0>=1)%>%
  ggplot(aes(
    y=hav_d,
    x=as.factor(N0),
    #fill=variable
  )) +
  geom_boxplot() +
  #geom_point(shape = 21, size = 3.1, fill = "grey") +
  #geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  #theme(legend.pos = "none") +
  ylab(expression(italic(d)["hav"])) + 
  xlab("Training samples within 1 Km") ->
  res.0
ggsave(res.0, file = "res.N0.pdf", w=4.5, h=4)