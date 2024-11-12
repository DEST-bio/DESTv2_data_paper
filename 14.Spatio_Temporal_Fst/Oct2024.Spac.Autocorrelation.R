#### Need geospatial R

#library(INLA); 
library(ggplot2); library(ggregplot)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ape)
library(nlme)
library(MuMIn)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")
samps %>%
  group_by(pop1=city) %>%
  summarize(lat.m = mean(lat),
            long.m = mean(long),
            continent=unique(continent)) -> lats

load("../ti.ob.1y.Rdata")
setDT(ti.ob.1y)
ti.ob.1y %<>%
  group_by(V1,V2,pop1,year1,year2) %>%
  slice_head()

left_join(ti.ob.1y, lats) -> ti.ob.1y

ti.ob.1y %>%
  filter(pop1 != "Providence") %>%
  filter(pop1 != "Yesiloz") %>%
  filter(year1 == 2015)->
  ti.ob.1y.clean
  
ti.ob.1y.clean %>%
  group_by(pop1, long.m, lat.m) %>%
  summarize(FST.m = mean(FST),
            T.mean.m = mean(T.mean, na.rm = T),
            T.max.m = mean(T.max, na.rm = T),
            T.min.m = mean(T.min, na.rm = T),
            T.var.m = mean(T.var, na.rm = T),
            Tn.below5.m = mean(Tn.below5, na.rm = T)
  ) ->
  ti.ob.means

ti.ob.means.DIST <- as.matrix(dist(cbind(ti.ob.means$long.m, 
                                      ti.ob.means$lat.m)))

ti.ob.means.DIST <- 1/ti.ob.means.DIST
diag(ti.ob.means.DIST) <- 0

#ti.ob.1y.DIST.inv[is.infinite(ti.ob.1y.DIST.inv)] <- 0

ti.ob.means.DIST[1:5, 1:5]

Moran.I(ti.ob.means$FST.m, ti.ob.means.DIST)
Moran.I(ti.ob.means$T.mean.m, ti.ob.means.DIST)
Moran.I(ti.ob.means$T.min.m, ti.ob.means.DIST)

Moran.I(ti.ob.means$T.max.m, ti.ob.means.DIST)
Moran.I(ti.ob.means$T.var.m, ti.ob.means.DIST, na.rm = TRUE)
Moran.I(ti.ob.means$Tn.below5.m, ti.ob.means.DIST)

####
base <- gls( FST.m ~ Tn.below5.m , 
                      data = ti.ob.means )

Gauss.autocor <- gls( FST.m ~ Tn.below5.m , 
                       correlation = corGaus(form = ~long.m + lat.m), 
                       data = ti.ob.means )
Exp.autocor <- gls( FST.m ~ T.var.m , 
                      correlation = corExp(form = ~long.m + lat.m), 
                      data = ti.ob.means )
Sph.autocor <- gls( FST.m ~ T.var.m , 
                      correlation = corSpher(form = ~long.m + lat.m), 
                      data = ti.ob.means )
Rat.autocor <- gls( FST.m ~ T.var.m , 
                      correlation = corRatio(form = ~long.m + lat.m), 
                      data = ti.ob.means )

model.sel(base, 
          Gauss.autocor, 
          Exp.autocor, 
          Sph.autocor, 
          Rat.autocor)

summary(base)
summary(Gauss.autocor)

###
base <- gls( FST.m ~ Tn.below5.m , 
             data = ti.ob.means )

Gauss.autocor <- gls( FST.m ~ Tn.below5.m , 
                      correlation = corGaus(form = ~long.m + lat.m), 
                      data = ti.ob.means )
Exp.autocor <- gls( FST.m ~ Tn.below5.m , 
                    correlation = corExp(form = ~long.m + lat.m), 
                    data = ti.ob.means )
Sph.autocor <- gls( FST.m ~ Tn.below5.m , 
                    correlation = corSpher(form = ~long.m + lat.m), 
                    data = ti.ob.means )
Rat.autocor <- gls( FST.m ~ Tn.below5.m , 
                    correlation = corRatio(form = ~long.m + lat.m), 
                    data = ti.ob.means )

model.sel(base, 
          Gauss.autocor, 
          Exp.autocor, 
          Sph.autocor, 
          Rat.autocor)

