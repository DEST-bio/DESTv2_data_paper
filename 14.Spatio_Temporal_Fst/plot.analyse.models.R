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
library(FactoMineR)
library(factoextra)
library(segmented)
library(gmodels)
library(tibble)
library(tidybulk)

### wheather
#nasapower <- get(load("/netfiles/nunezlab/Drosophila_resources/NASA_power_weather/DEST2.0/nasaPower.allpops.Rdata"))

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

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
save(ti.ob.1y, file = "ti.ob.1y.Rdata")
####
#load("fst.winter.1y.Rdata")
####

setDT(ti.ob.1y)
ti.ob.1y %<>%
  filter(pop1 != "Providence")
###
unique(ti.ob.1y$pop1)
unique(ti.ob.1y$pop2)

### broken stick regression
samps %>%
  group_by(pop1=city) %>%
  summarize(lat.m = mean(lat)) -> lats

left_join(ti.ob.1y, lats) -> ti.ob.1y

ti.ob.1y %>%
group_by(pop1) %>%
summarize(
mT = mean(T.mean , na.rm = T),
meanL = mean(lat.m),
meanF = ci(FST)[1],
lF = ci(FST)[2],
hF = ci(FST)[3],
tmax = mean(T.max)
) -> summaries_lst_fst 

my.lm <- lm(meanF ~ meanL, data = summaries_lst_fst)
summary(my.lm)

my.seg <- segmented(my.lm, 
                    seg.Z = ~ meanL, 
                    psi = 50.33469)

summary(my.seg)
my.seg$psi

my.fitted <- fitted(my.seg)
my.model <- data.frame(meanL = summaries_lst_fst$meanL, 
						meanF = my.fitted)


#### plot here!
ggplot()+
geom_vline(xintercept = 50.335, linetype = "dashed") +
geom_line(data = my.model, 
aes(x=meanL,y=(meanF)), 
colour = "tomato", linewidth = 1.9) +
geom_errorbar(width = 0.1, data = summaries_lst_fst,
aes(
x=meanL,
ymin=(lF),
ymax=(hF),
)) +
geom_point(shape = 21, 
size = 4, 
data = summaries_lst_fst,
aes(
x=meanL,
y=((meanF)),
fill = mT
)) + 
  theme_bw() + ylim((0.0), (0.1)) +
  scale_fill_gradient2(low= "steelblue", 
                       high="firebrick4", 
                       midpoint = 15) ->
datawbkstick
ggsave(datawbkstick, file = "datawbkstick.pdf", h = 3, w= 4)



#### Extreme analysis
ti.ob.1y %>%
ggplot(
aes(
x=lat.m,
y=(FST),
group=pop1
)) + geom_boxplot(width = 0.25, outlier.colour="red", outlier.shape=8) +
geom_vline(xintercept = 50.335, linetype = "dashed") +
theme_bw()  ->
FSTplo
ggsave(FSTplo, file = "FSTplo.pdf", h = 4, w= 5)

####

####
cor.test(logit(abs(ti.ob.1y$FST)), ti.ob.1y$day_diff )
cor.test(logit(abs(ti.ob.1y$FST)), ti.ob.1y$lat.m )

cor.test(ti.ob.1y$T.mean, ti.ob.1y$lat.m )

ti.ob.1y %>%
  filter(pop1 != "Yesiloz") -> noTur
cor.test(logit(abs(noTur$FST)), noTur$lat.m )
cor.test(logit(abs(noTur$FST)), noTur$T.mean )


ti.ob.1y %>%
  filter(pop1 == "Yesiloz") -> Tur
cor.test(logit(abs(Tur$FST)), Tur$lat.m )
cor.test(logit(abs(Tur$FST)), Tur$T.mean )

model <- aov(logit(abs(ti.ob.1y$FST)) ~ lat.m, data = ti.ob.1y)
model$coefficients
summary(model)
Anova(model)


### Make plot of FST and latitude
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

### Permute findings
perms_overw = 
foreach(i=1:500, .combine = "rbind")%do%{
    data.frame(
      i = i,
      cor=cor.test(ti.ob.1y$FST, sample(ti.ob.1y$lat.m) )$estimate 
    )
}

rbind(
data.frame(perm.stat = "all", i = 0, cor= cor.test(~FST+lat.m, data = ti.ob.1y)$estimate),
data.frame(perm.stat = "less50", i = 0, cor= cor.test(~FST+lat.m, data = filter(ti.ob.1y, lat.m < 50.3))$estimate),
data.frame(perm.stat = "more50", i = 0, cor= cor.test(~FST+lat.m, data = filter(ti.ob.1y, lat.m > 50.3))$estimate),
data.frame(perm.stat = "perm", perms_overw)) -> perm.test
## plot permutations
ggplot() +
  geom_density(data = filter(perm.test, perm.stat == "perm"),
               aes(cor), fill = "grey"
               ) + 
  geom_vline(data = filter(perm.test, perm.stat == "all"),
             aes(xintercept=cor), color = "purple") +
    geom_vline(data = filter(perm.test, perm.stat == "less50"),
             aes(xintercept=cor), color = "blue") +
      geom_vline(data = filter(perm.test, perm.stat == "more50"),
             aes(xintercept=cor), color = "red") +           
  theme_bw() ->
  perm.plots
ggsave(perm.plots, file = "perm.plots.pdf", w = 4, h = 4)




####
ti.ob.1y %>%
  filter(pop1 == "Yesiloz") %>%
  mutate(year_l = paste(year1,year2)) %>%
  ggplot(aes(
    x=T.mean,
    y=(abs(.$FST)),
    fill=as.numeric(T.mean)) 
  ) +
  geom_point(size = 3, shape = 21) +
  ylim(0, 0.1) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_wrap(~year_l=="2020 2021", scale = "free_x") +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot.temp
ggsave(lat.fst.plot.temp, 
       file = "Yesiloz.lat.fst.plot.temp.pdf", 
       w = 5.5, h = 3.5)

#### Eco PCA variables
ti.ob.1y %>%
dplyr::select(T.mean,T.var,T.min,T.max,Tn.below5,Tn.above32) %>%
PCA(graph = FALSE) -> eco.PCA

fviz_pca_var(eco.PCA, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
             ) -> eco.PCA.vars
  ggsave(eco.PCA.vars, 
       file = "eco.PCA.vars.pdf", 
       w = 6, h = 3)
          

eco.PCA$ind$coord %>%
as.data.frame %>%
mutate(ti.ob.1y) -> eco.fst.dat

eco.fst.dat %>% 
  dplyr::select(Dim.1, T.min, FST, Tn.below5) %>%
  reshape2::melt(id = c("Dim.1","FST") ) %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
    fill =FST
  )) + geom_point(alpha = 0.1, shape = 21, size = 3) +
  scale_fill_gradient2(mid = "white", high = "springgreen", midpoint = 0) +
  facet_wrap(~variable) ->
  dim1.cold
ggsave(dim1.cold, file ="dim1.cold.pdf")

eco.fst.dat %>% 
  dplyr::select(Dim.1, T.min, FST, Tn.below5) %>%
  reshape2::melt(id = c("Dim.1","FST") ) %>%
  ggplot(aes(
    x=FST ,
    y=value,
    fill =Dim.1
  )) + geom_point(alpha = 0.1, shape = 21, size = 3) +
  scale_fill_gradient2(mid = "white", high = "springgreen", midpoint = 0) +
  facet_wrap(~variable) ->
  dim1.cold.fst
ggsave(dim1.cold.fst, file ="dim1.cold.fst.pdf")


#eco.fst.dat %>% 

  filter(eco.fst.dat, pop1 != "Yesiloz") %>%
    filter(Tn.below5 > 1) %>%
  dplyr::select(Dim.2, T.max, Tn.below5, FST) %>%
  reshape2::melt(id = c("Dim.2","FST","Tn.below5") ) %>%
  ggplot(aes(
    x=Tn.below5 ,
    y=value,
    fill =FST
  )) + geom_point(alpha = 0.1, shape = 21, size = 3) +
    geom_smooth(method = "lm") +
  scale_fill_gradient2(mid = "white", high = "springgreen", midpoint = 0) +
  facet_wrap(~variable) ->
  dim2.hot.fst
ggsave(dim2.hot.fst, file ="dim2.hot.fst.pdf")


### Examine variable associations
#cor.test(~FST+T.min, data = eco.fst.dat)
cor.test(~FST+T.min, data = filter(eco.fst.dat, pop1 != "Yesiloz"))
## T.max!
cor.test(~FST+T.max, data = filter(eco.fst.dat, pop1 != "Yesiloz"))

#### --> Tn.below5
cor.test(~FST+Tn.below5 , data = filter(eco.fst.dat, pop1 != "Yesiloz"))
eco.fst.dat %>%
  filter(Tn.below5 > 1) %>%
  filter(pop1 != "Yesiloz") -> tnb5
cor.test(~FST+Tn.below5 , data = tnb5)
cor.test(~FST+Tn.below5 , data = filter(tnb5, lat.m > 50))
cor.test(~FST+Tn.below5 , data = filter(tnb5, lat.m < 50))

tnb5 %>%
ggplot(aes(
x=Tn.below5,
y=((FST)),
fill =T.min,
)) + geom_point(size = 3, aes(shape=lat.m > 50.3)) +
geom_smooth(method = "lm", color = 'black', linetype = "solid") +
theme_bw() + ylim(0,0.1)+
  scale_fill_gradient2(low="darkblue", high="steelblue", midpoint = -10) +
  scale_shape_manual(values = 21:22) ->
eco.fst 
ggsave(eco.fst, file = "eco.fst.pdf", w = 6, h = 3)


####
ti.ob.1y %>%
  filter(pop1 != "Yesiloz") %>%
  ggplot(aes(
    x=Tn.below5,
    y=logit(abs(FST)),
    fill=as.numeric(T.mean)) 
  ) +
  geom_point(size = 3, shape = 21) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(~lat.m > 50) +
  scale_fill_gradient2(low="steelblue", high="firebrick4", midpoint = 15) ->
  lat.fst.plot.tempvar
ggsave(lat.fst.plot.tempvar, 
       file = "lat.fst.plot.tempvar.pdf", 
       w = 6, h = 3)


##### Map
samps %>%
  group_by(pop1=city) %>%
  summarize(lat.m = mean(lat),
            long.m = mean(long)) -> lats

left_join(ti.ob.1y, lats) -> ti.ob.1y

ti.ob.1y %>%
  group_by(pop1, lat.m, long.m) %>%
  summarise(FST_m = mean(FST)) -> FST_sum

world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1
  ) + theme_classic() -> base_world


base_world + 
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-140, 41.00), ylim = c(20.00, 68.00), expand = FALSE) + 
  geom_point(
    data = FST_sum,
    aes(x=long.m,
        y=lat.m,
        fill = FST_m), size = 3.5, shape = 21
  ) + scale_fill_gradient2(mid = "white", high = "springgreen", midpoint = 0) ->
  map.overfst
ggsave(map.overfst, 
       file = "map.overfst.pdf", 
       w = 9, h = 3)

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


