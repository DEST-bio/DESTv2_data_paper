
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
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(forcats)
library(FactoMineR)

#### samps
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv"

samps <- fread(meta_git)
setDT(samps)

##### ---> Part 1: DEST 1.0 Model
###
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

cor.test(~Region_havd+City_havd, data =DEST1.comps)

###
###
load("DEST2.0.model.predictions.Rdata")

rbind(
mutate(hav.dist.obj.top.pred[,c("continent","hav_d")], model = "DEST 1.0"), 
mutate(
dplyr::select(gims.out.v2, continent=cont,hav_d )[,-1], model = "DEST 2.0")) %>%
group_by(continent, model) %>%
summarize(mean.d = ci(hav_d)[1],
		  lci= ci(hav_d)[2],
		  uci= ci(hav_d)[3],
		  m.sd = sd(hav_d)
) -> mod.summaries

d1.med = median(filter(mod.summaries, model == "DEST 1.0")$mean.d)
d2.med = median(filter(mod.summaries, model == "DEST 2.0")$mean.d)

mod.summaries %>%
  ggplot(aes(
    x=continent,
    y=mean.d,
    ymin = mean.d - m.sd,
    ymax = mean.d + m.sd,
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
  coord_flip() ->  DAPC.v1.v2
  #scale_y_continuous(trans='log10') 

ggsave(DAPC.v1.v2, file = "DAPC.v1.v2.model.performance.pdf", w= 5, h = 4)

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



