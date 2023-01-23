library(foreach)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)
library(readxl)


### Load metadata
#setwd("/Users/alanbergland/Documents/GitHub")
samps <- fread("../DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")

### only use samps collected in a small period
samps <- samps
table(samps$set)

### only use localities with more than one sample per year
samps.ag <- samps[,list(.N, daysDelta=max(jday) - min(jday)), list(locality, year)]

setkey(samps, locality, year)
setkey(samps.ag, locality, year)
samps <- merge(samps, samps.ag)

table(samps$set)

####
### which samples were in the core20
drosrtec <- as.data.table(read_excel("../DESTv2/populationInfo/OriginalMetadata/elife-67577-supp1-v2.xlsx"))
drosrtec[,Sample:=gsub("PA_sc", "PA_st", Sample)]

samps <- merge(samps, drosrtec[,c("Sample", "Core20")], by.x="sampleId_orig", by.y="Sample", all.x=T)
samps[Core20=="yes", Core20:=T]
samps[is.na(Core20) | Core20=="no", Core20:=F]

table(samps$Core20, samps$set)
table(samps[N>1]$Core20, samps[N>1]$set)
table(samps[N>1 & daysDelta>45]$Core20, samps[N>1 & daysDelta>45]$set)
samps.ag[N>1 & daysDelta>45] ->
  initial_samp_filter

##########
### Load NASA power
load("/project/berglandlab/alan/nasaPower.allpops.Rdata")

weather.explo = 
foreach(i=1:dim(initial_samp_filter)[1], .combine = "rbind")%do%{

  message(paste(i, initial_samp_filter$locality[i]))
################
samps %>%
  filter(locality == initial_samp_filter$locality[i]) %>%
  filter(year == initial_samp_filter$year[i]) ->
  tmp.samps

foreach(j=1:dim(tmp.samps)[1], .combine = "rbind")%do%{
  
focal.samp = tmp.samps$sampleId[j]
core20.stat = tmp.samps$Core20[j]
message(paste(i, initial_samp_filter$locality[i], focal.samp))

##############

month(as.Date(strsplit(focal.samp, "_")[[1]][5], format = "%Y-%m-%d")) -> j.tmp.month

if(j.tmp.month %in% 5:7){season = "spring"} else
  if(j.tmp.month %in% 9:11){season = "fall"}
message(season)

if(season == "spring"){time.adjustment = +14} else
  if(season == "fall"){time.adjustment = -14}
message(time.adjustment)



power.dt %>%
  filter(sampleId == focal.samp) %>% 
  separate(date, into = c("Date.power", "hour.power"), sep = "\ ") %>% 
  mutate(Date.power = as.Date(Date.power, format = "%Y-%m-%d")) %>%
  group_by(sampleId, Date.power) %>% 
  summarize(min.T = min(T2M),
            max.T = max(T2M),            
            mean.T = mean(T2M),
            var.T = var(T2M)
            ) %>% 
  separate(sampleId, remove = F, into = c("cont","stat", "city", "rep", "date.samp"), sep = "_") %>%
  mutate(date.samp = as.Date(date.samp, format = "%Y-%m-%d")) %>%
  mutate(date_boundary.winter =  as.Date(date.samp)+time.adjustment) %>%
  mutate(date_boundary.summer =  as.Date(date.samp)+((time.adjustment)*-1)) ->
  dat.power

### Get boundary conditions  
if(season == "spring"){
  dat.power %>%
    filter(Date.power <= date_boundary.winter & Date.power >=  date.samp) ->
    time.frame.condit.winter
  dat.power %>%
    filter(Date.power >= date_boundary.summer & Date.power <=  date.samp) ->
    time.frame.condit.summer
} else
  if(season == "fall"){
    dat.power %>%
      filter(Date.power >= date_boundary.winter & Date.power <=  date.samp) ->
      time.frame.condit.winter
    dat.power %>%
      filter(Date.power <= date_boundary.summer & Date.power >=  date.samp) ->
      time.frame.condit.summer
  }

### Frost test
time.frame.condit.winter %>%
  summarize(frost.test = sum(min.T <= 0)) -> freeze.test

### 30 pick by WHO --> https://www.who.int/india/heat-waves
### Heat wave test
time.frame.condit.summer %>%
  summarize(heatwave.test = sum(max.T >= 35)) -> heat.test

### Get 15d temperatures
### This collect environemtal data fifteen days prior to collection regarldess of seasonal "bound"
dat.power %>%
  filter(Date.power <= date.samp & Date.power >= date.samp-14) ->
  temp.15d
##########

paste(sep ="_",
strsplit(focal.samp, "_")[[1]][1],
strsplit(focal.samp, "_")[[1]][2],
strsplit(focal.samp, "_")[[1]][3]
) -> locale.tmp

data.frame(
  sampleId=focal.samp,
  locality=locale.tmp,
  season=season,
  time.adjustment=time.adjustment,
  freeze.test=freeze.test$frost.test,
  heat.test=heat.test$heatwave.test
) %>%
  mutate(Tmean.15 = mean(temp.15d$mean.T),
         Tvar.15 = var(temp.15d$mean.T),
         Tmax.15 = max(temp.15d$mean.T),
         Tmin.15 = min(temp.15d$mean.T),
         PropTmin5.15 = sum(temp.15d$min.T < 5),
         PropTmax32.15 = sum(temp.15d$max.T > 32),
         Core20_sat = core20.stat
         )-> out.tmp

return(out.tmp)

} ### k loop
} ### i loop


########
### apply frost filter
weather.explo %>%
  filter(Core20_sat == FALSE) %>%
  group_by(sampleId) %>%
  mutate(year = year(as.Date(strsplit(sampleId, "_")[[1]][5]))) %>%  
  filter(freeze.test == 0 & heat.test <= 3) %>%
  .$sampleId ->
  pass.filters

weather.explo %>%
  filter(sampleId %in% pass.filters) %>% 
  group_by(sampleId) %>%
  mutate(  year = year(strsplit(sampleId, "_")[[1]][5]) ) %>%
  group_by(locality, year, loc.y = paste(locality, year, sep = "_")) %>% 
  summarize(N = n()) %>%
  filter(N >= 2) -> 
  filtered_locs

weather.explo %>%
  group_by(sampleId) %>%
  mutate(year = year(strsplit(sampleId, "_")[[1]][5]) ) %>%
  group_by(locality, year, loc.y = paste(locality, year, sep = "_")) %>% 
  filter(sampleId %in% pass.filters) %>%
  filter(loc.y %in% filtered_locs$loc.y) ->
  seasonal_pairs.cands

weather.explo %>%
  filter(Core20_sat == TRUE ) %>% 
  group_by(sampleId) %>%
  mutate(year = year(strsplit(sampleId, "_")[[1]][5]) ) %>%
  group_by(locality, year, loc.y = paste(locality, year, sep = "_")) ->
  samps.core20

###
seasonal_pairs.cands %>%
  filter(season == "spring") %>%
  group_by(loc.y) %>%
  slice_head(n=1) -> springs

seasonal_pairs.cands %>%
  filter(season == "fall") %>%
  group_by(loc.y) %>%
  slice_tail(n=1) -> falls

###
rbind(springs, falls) %>%  
  dcast(loc.y~season, value.var = "sampleId") %>%
  .[complete.cases(.),] ->
  final.filter

final.samps = c(final.filter$fall, final.filter$spring)

rbind(springs, falls) %>%  
  filter(sampleId %in% final.samps) ->
  DEST2.seasonals

######
rbind(DEST2.seasonals, samps.core20) ->
  DEST2.seasonals.plusCore20
  
DEST2.seasonals.plusCore20 %>%
  group_by(Core20_sat) %>% summarize(N=n())

### Flipped samples
DEST2.seasonals.plusCore20 %>%
  left_join(dplyr::select(samps, sampleId, sampleId_orig), by ="sampleId") %>%
  mutate(flipped = case_when(
  sampleId_orig %in% c(paste(c("KA_to_14", "MA_la_12", "MI_bh_14", "CA_es_12"), "_fall", sep = ""),
                       paste(c("KA_to_14", "MA_la_12", "MI_bh_14", "CA_es_12"), "_spring", sep = "")
  ) ~ "fliped",
  TRUE ~ "not.flipped"
)) -> DEST2.seasonals.plusCore20.flip

DEST2.seasonals.plusCore20.flip %>% 
  filter(loc.y == "US_Wis_Cro_2012")

DEST2.seasonals.plusCore20.flip$flipped %>% table

####
DEST2.seasonals.plusCore20.flip %>% 
  group_by(sampleId) %>%
  mutate(month = month(strsplit(sampleId, "_")[[1]][5]) ) %>%
  arrange(month) %>%
  group_by(loc.y) %>%
  mutate(order= 1:2) %>% as.data.frame() %>% 
  dcast(loc.y~order, value.var = "Tmean.15") %>%
  mutate(T.dir = `2`-`1`) %>%
  mutate(delta.T.sign = sign(T.dir)) %>%
  mutate(delta.T.mag = case_when( abs(T.dir) >= 0 & abs(T.dir) <= 1 ~ "Flat",
                                  abs(T.dir) > 1 &  abs(T.dir) <= 3 ~ "Shallow",
                                  abs(T.dir) >= 3 ~ "Steep")) ->
  seas.meta

#### final merge
DEST2.seasonals.plusCore20.flip %>%
  group_by(sampleId) %>%
  mutate(date.samp = (strsplit(sampleId, "_")[[1]][5]) ) %>%
  left_join(seas.meta) ->
  DEST2.seasonals.plusCore20.flip.met

####
save(DEST2.seasonals.plusCore20.flip.met, file = "DEST2.seasonals.plusCore20.flip.met.Rdata")
load("DEST2.seasonals.plusCore20.flip.met.Rdata")

### plot
DEST2.seasonals.plusCore20.flip.met %>%
  ggplot(aes(
    x=month(date.samp),
    y=Tmean.15,
    fill = season,
    group = loc.y,
    shape = flipped
  )) +
  geom_line(color = "black", linetype = "dashed") +
  geom_point(size = 2) +
  scale_fill_manual(values = c("red","blue")) +
  scale_shape_manual(values = 21:22) +
  theme_bw() +
  #geom_hline(yintercept = 3) +
  facet_grid(Core20_sat~delta.T.mag*sign(T.dir)) ->
  seas.pairs

ggsave(seas.pairs, file = "seas.pairs.pdf", h =4, w = 8)

#### QC conditions
qc.dat <- get(load("/scratch/yey2sn/seasonlaity_choice/sample_Assesment/DEST2.o.all.dmel.delim.Rdata"))
sim.contam <- get(load("./sim.contam.joint.Rdata"))

names(qc.dat)[3] = "sampleId"

left_join(DEST2.seasonals.plusCore20.flip.met, qc.dat, by = "sampleId") %>% 
  left_join(dplyr::select(samps, sampleId, nFlies), by = "sampleId") %>%
  filter(bam_slice == "Dmel") %>%
  filter(metric == "coverage_across_reference") %>%
  group_by(sampleId, nFlies) %>% 
    #summarize(m.var = mean(var)) %>% 
  group_by(sampleId) %>%
  mutate(top = (var*(2*nFlies)-1),
         bottom = (var+2*nFlies)) %>%
  mutate(Neff = top/bottom) -> Neff.dat

Neff.dat %>% 
  group_by(sampleId) %>%
  summarize(me.neff = median(Neff)) %>%
  arrange(me.neff)

# MA_Tan_Lar_1_2021-09-06   0.948
#2 FR_Hau_Civ_1_2018-09-15   1.02 
#3 UA_Ode_Ode_1_2019-06-15  11.1  

Neff.dat %>%
  ggplot(aes(
    x=sampleId,
    y=Neff,
    fill=season
  )) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 20) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) ->
  neff.box

ggsave(neff.box, file = "neff.box.pdf", w = 9 , h =2.5)
####
left_join(DEST2.seasonals.plusCore20.flip.met, sim.contam, by = "sampleId") ->
  contam

contam %>% 
  group_by(sampleId) %>%
  summarize(m.sim.c = median(m.sim)) %>%
  arrange(-m.sim.c) %>%
  filter(m.sim.c > 0.05)

contam %>%
  ggplot(aes(
    x=sampleId,
    y=m.sim,
    fill=season
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) ->
  sim.contam.box
ggsave(sim.contam.box, file = "sim.contam.box.pdf", w = 9 , h =2.5)


#### tables
contam %>% 
  group_by(sampleId) %>%
  summarize(m.sim.c = median(m.sim)) %>%
  arrange(-m.sim.c)  -> c.samps

Neff.dat %>% 
  group_by(sampleId) %>%
  summarize(me.neff = median(Neff)) %>%
  arrange(me.neff) -> neff.samps

full_join(c.samps, neff.samps) %>% 
  melt(id = "sampleId") %>% 
  ggplot(aes(
    x=as.character(sampleId),
    y=value
  )) +
  geom_boxplot() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  facet_wrap(~variable, ncol =1, scales = "free_y") ->
  qc.fig.dmel

ggsave(qc.fig.dmel, file = "qc.fig.dmel.pdf", w = 8, h =4)

full_join(c.samps, neff.samps) %>%
  group_by(sampleId) %>%
  filter(me.neff < 20 | m.sim.c > 0.05) %>%
  dplyr::select(sampleId, me.neff, m.sim.c)
