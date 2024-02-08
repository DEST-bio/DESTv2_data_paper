### Plot FST analysis spac. temp

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

### wheather
#nasapower <- get(load("/netfiles/nunezlab/Drosophila_resources/NASA_power_weather/DEST2.0/nasaPower.allpops.Rdata"))

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

#### Many winters
#### Many winters
#### Many winters
#### Many winters
#### Many winters

load("../fst.winter.newTimePops.Rdata")
setDT(fst.winter.between) 

fst.winter.between %<>%
  filter(pop1 != "Karensminde")


####
#anova(lm(logit(abs(FST)) ~ pop1*day_diff+T.min+T.max+T.mean, data = fst.winter.between))
#lm(logit(abs(FST)) ~ pop1*day_diff, data = fst.winter.between) -> mod.tim.loc

cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.max)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.min)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.mean)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$day_diff)

foreach(i=unique(fst.winter.between$pop1),
        .combine = "rbind")%do%{
          
          message(i)
          tmp = fst.winter.between %>%
            filter(pop1 == i)
    
      cor.test(logit(abs(tmp$FST)), tmp$T.mean)-> tp
      cor.test(logit(abs(tmp$FST)), tmp$T.min)-> tmin
      cor.test(logit(abs(tmp$FST)), tmp$T.max)-> tmax
      cor.test(logit(abs(tmp$FST)), tmp$day_diff)-> timep
      
      data.frame(
        pop=i,
        meanT.p = tp$p.value,
        minT.corr = tmin$estimate,
        minT.p = tmin$p.value,
        maxT.corr = tmax$estimate,
        maxT.p = tmax$p.value,
        meanT.corr = tp$estimate,
        Time.p = timep$p.value,
        Time.corr = timep$estimate
      )
      
        }

####
fst.winter.between %>%
  dplyr::select(#residual.t.loc,
    #ydelta.test,
    day_diff,
    pop1,
    T.mean,
    FST,
    #T.min,
    #T.max
  ) %>%
  melt(id=c("FST","pop1"
  ) ) %>% 
  ggplot(
    aes(
      x=(value),
      y=logit(abs(FST)),
      color = pop1,
    )
  ) +
  ylab(expression(paste(#"residual", 
    F[ST]))) +
  xlab("value between 2 samples") +
  geom_point(alpha = 0.1) +
  facet_grid(.~variable,
             #ncol = 2,
             scales = "free_x"
  ) +
  theme_bw() +
  geom_smooth(method = "lm", se =T, color = "black")->
  Twinbet.ti.plot
ggsave(Twinbet.ti.plot, 
       file = "Twinbet.ti.plot.pdf", 
       w = 7, h = 3)

########
########
########

fst.winter.between %>%
  dplyr::select(#residual.t.loc,
    #ydelta.test,
    day_diff,
    pop1,
    T.mean,
    FST,
    #T.min,
    #T.max
  ) %>%
  melt(id=c("FST","pop1"
  ) ) %>% 
  ggplot(
    aes(
      x=(value),
      y=logit(abs(FST)),
      color = pop1,
    )
  ) +
  ylab(expression(paste(#"residual", 
    F[ST]))) +
  xlab("value between 2 samples") +
  #geom_point(alpha = 0.1) +
  facet_wrap(.~variable,
             ncol = 2,
             scales = "free_x"
  ) +
  theme_bw() +
  geom_smooth(method = "lm", se =T)->
  Twinbet.ti.plot.pop
ggsave(Twinbet.ti.plot.pop, 
       file = "Twinbet.ti.plot.pop.pdf", 
       w = 4, h = 3)

########
########
######## ---> one year analysis
######## ---> one year analysis
######## ---> one year analysis
######## ---> one year analysis
######## ---> one year analysis
######## ---> one year analysis


###
fld.1y = "/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/weather_fst_1y"

###
ti.f.1y = system(paste("ls ", fld.1y, sep = "" ), intern = T)

##
ti.ob.1y =
  foreach(i=ti.f.1y,
          .combine = "rbind"
  )%do%{
    tmp = get(load(paste(fld.1y,i, sep = "/")))
    return(tmp)
  }

#load("fst.winter.1y.Rdata")
setDT(ti.ob.1y)
ti.ob.1y %<>%
  filter(pop1 != "Providence")
###
unique(ti.ob.1y$pop1)
unique(ti.ob.1y$pop2)

###
samps %>%
  group_by(pop1=city) %>%
  summarize(lat.m = mean(lat)) -> lats

left_join(ti.ob.1y, lats) -> ti.ob.1y
####
cor.test(logit(abs(ti.ob.1y$FST)), ti.ob.1y$day_diff )
cor.test(logit(abs(ti.ob.1y$FST)), ti.ob.1y$lat.m )

ti.ob.1y %>%
  filter(pop1 != "Yesiloz") -> noTur
cor.test(logit(abs(noTur$FST)), noTur$lat.m )
cor.test(logit(abs(noTur$FST)), noTur$T.mean )


ti.ob.1y %>%
  filter(pop1 == "Yesiloz") -> Tur
cor.test(logit(abs(Tur$FST)), Tur$lat.m )
cor.test(logit(abs(Tur$FST)), Tur$T.mean )

model <- aov(logit(abs(ti.ob.1y$FST)) ~ lat.m, data = ti.ob.1y)
summary(model)
Anova(model)

###
###
ti.ob.1y %>%
  ggplot(aes(
    x=lat.m,
    y=logit(abs(ti.ob.1y$FST)),
    fill=T.mean,
    shape = pop1=="Yesiloz") 
  ) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  scale_shape_manual(values = c(21,22)) +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot
ggsave(lat.fst.plot, 
       file = "lat.fst.plot.pdf", 
       w = 6, h = 5)

ti.ob.1y %>%
  ggplot(aes(
    x=T.mean,
    y=logit(abs(ti.ob.1y$FST)),
    fill=T.mean,
    shape = pop1=="Yesiloz") 
  ) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(pop1=="Yesiloz"~.) +
  scale_shape_manual(values = c(21,22)) +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot.temp
ggsave(lat.fst.plot.temp, 
       file = "lat.fst.plot.temp.pdf", 
       w = 6, h = 4)


ti.ob.1y %>%
  ggplot(aes(
    x=day_diff,
    y=logit(abs(ti.ob.1y$FST)),
    fill=T.mean,
    shape = pop1=="Yesiloz") 
  ) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(pop1=="Yesiloz"~.) +
  scale_shape_manual(values = c(21,22)) +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot.day_diff
ggsave(lat.fst.plot.day_diff, 
       file = "lat.fst.plot.day_diff.pdf", 
       w = 6, h = 4)

### Comparing variance in temperature and in FST

ti.ob.1y %>%
  group_by(pop1) %>%
  summarize(mean.lat = mean(lat.m), varT = var(T.mean, na.rm = T)) %>%
  ggplot(aes(
    y=varT,
    x=mean.lat
  )) +
  geom_point() ->
  var.plot

ggsave(var.plot, 
       file = "var.plot.pdf", 
       w = 6, h = 4)


