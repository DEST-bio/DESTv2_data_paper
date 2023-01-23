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

samps %>%
  group_by(sampleId) %>%
  mutate(date = as.Date(strsplit(sampleId, "_")[[1]][5], format = "%Y-%m-%d")) %>% 
  mutate(month = month(date)  ) %>% 
  filter(month %in% 9:11 ) ->
  fall.samps


fall.samps %>%
  mutate(loc.y = paste(locality, year, sep = "_")) %>%
  group_by(loc.y) %>%
  summarize(N = n()) %>%
  filter(N >= 2) ->
  winter.loc.y

fall.samps %>% 
  mutate(loc.y = paste(locality, year, sep = "_")) %>%
  filter(loc.y %in% winter.loc.y$loc.y) ->
  fall.samps.filtered

##########
### Load NASA power
load("/project/berglandlab/alan/nasaPower.allpops.Rdata")

weather.freeze = 
foreach(i=1:dim(fall.samps.filtered)[1], .combine = "rbind")%do%{

  focal.samp = fall.samps.filtered$sampleId[i]
  
##############

year(as.Date(strsplit(focal.samp, "_")[[1]][5], format = "%Y-%m-%d")) -> j.tmp.year
month(as.Date(strsplit(focal.samp, "_")[[1]][5], format = "%Y-%m-%d")) -> j.tmp.month

power.dt %>%
  filter(sampleId == focal.samp) %>% 
  filter(YEAR == j.tmp.year ) %>%
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
  mutate(date_boundary.before =  as.Date(date.samp)-14) %>%
  mutate(date_boundary.after =  as.Date(date.samp)+14) ->
  dat.power

###
dat.power %>%
  filter(Date.power >=  date.samp & Date.power <= date_boundary.after ) ->
  dat.after
dat.power %>%
  filter(Date.power <=  date.samp & Date.power >= date_boundary.before ) ->
  dat.before
###

### Get boundary conditions  
### Frost test
dat.after %>%
  summarize(frost.test.after = sum(min.T <= 0)) -> freeze.test.after
dat.before %>%
  summarize(frost.test.before = sum(min.T <= 0)) -> freeze.test.before

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
  season="frost.set",
  freeze.after=freeze.test.after$frost.test.after,
  freeze.before=freeze.test.before$frost.test.before
) %>%
  mutate(Tmean.15 = mean(temp.15d$mean.T),
         Tvar.15 = var(temp.15d$mean.T),
         Tmax.15 = max(temp.15d$max.T),
         Tmin.15 = min(temp.15d$min.T),
         PropTmin5.15 = sum(temp.15d$min.T < 5),
         PropTmax32.15 = sum(temp.15d$max.T > 32),
         )-> out.tmp

return(out.tmp)

} ### i loop

weather.freeze %>% 
  group_by(sampleId) %>%
  mutate(date = (as.Date(strsplit(sampleId, "_")[[1]][5], format = "%Y-%m-%d"))) %>%
  mutate(year = year(date)) %>%
  mutate(loc.y = paste(locality, year, sep = "_")) %>%
  filter(freeze.before < 1) %>%
  group_by(loc.y) %>%
  mutate(status = "pre.frost") %>%
  slice_tail(n=1) ->
  no.freeze

weather.freeze %>% 
  group_by(sampleId) %>%
  mutate(date = (as.Date(strsplit(sampleId, "_")[[1]][5], format = "%Y-%m-%d"))) %>%
  group_by(sampleId) %>%
  mutate(year = year(date)) %>%
  mutate(loc.y = paste(locality, year, sep = "_")) %>%
  filter(freeze.before >= 1) %>%
  group_by(loc.y) %>%
  mutate(status = "post.frost") %>%
  slice_head(n=1) ->
  post.freeze

rbind(no.freeze, post.freeze) %>%
  dcast(loc.y~status, value.var = "sampleId") %>% 
  .[complete.cases(.),] ->
  pre.post.frost.samps

###
weather.freeze %>%
  filter(sampleId %in% c(pre.post.frost.samps$post.frost,pre.post.frost.samps$pre.frost) ) ->
  samps.frost.final.set
####
save(samps.frost.final.set, file = "frost.final.set.Rdata")

load("./frost.final.set.Rdata")

samps.frost.final.set %>%
  group_by(sampleId) %>%
  mutate(city=strsplit(sampleId, "_")[[1]][2]) %>% 
  mutate(date = (as.Date(strsplit(sampleId, "_")[[1]][5], format = "%Y-%m-%d"))) %>%
  mutate(year = year(date)) %>%
  mutate(loc.y = paste(locality, year, sep = "_")) %>%
  ggplot(aes(
    x=(date),
    y=Tmin.15,
    fill = city,
    group = loc.y,
  )) +
  geom_hline(yintercept = 0) +
  geom_line(color = "black", linetype = "dashed") +
  geom_point(size = 2, shape =21) +
  #scale_fill_manual(values = c("red","blue")) +
  #scale_shape_manual(values = 21:22) +
  facet_wrap(~year, scales = "free_x") +
  theme_bw() ->
  freeze.pairs

ggsave(freeze.pairs, file = "freeze.pairs.pdf", h =4, w = 8)

####
qc.dat <- get(load("/scratch/yey2sn/seasonlaity_choice/sample_Assesment/DEST2.o.all.dmel.delim.Rdata"))
names(qc.dat)[3] = "sampleId"

sim.contam <- get(load("./sim.contam.joint.Rdata"))

left_join(samps.frost.final.set, qc.dat, by = "sampleId") %>% 
  left_join(dplyr::select(samps, sampleId, nFlies), by = "sampleId") %>%
  left_join(sim.contam) %>% 
  filter(bam_slice == "Dmel") %>%
  filter(metric == "coverage_across_reference") %>%
  group_by(sampleId, nFlies, m.sim) %>%  
  summarize(m.var = mean(var)) %>% 
  mutate(top = (m.var*(2*nFlies)-1),
         bottom = (m.var+2*nFlies)) %>%
  mutate(Neff = top/bottom) -> Neff.dat
