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

#######
args = commandArgs(trailingOnly=TRUE)
k= as.numeric(args[1]) 

### wheather
#nasapower <- get(load("/netfiles/nunezlab/Drosophila_resources/NASA_power_weather/DEST2.0/nasaPower.allpops.Rdata"))

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

###
fld = "/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/outfiles_1y"

###
ti.f = system(paste("ls ", fld, sep = "" ), intern = T)

nasa.power = get(load("Nasapower.winter.range.DEST.Rdata"))

####
#ti.ob =
#foreach(i=ti.f,
#.combine = "rbind"
#)%do%{
#tmp = get(load(paste(fld,i, sep = "/")))
#return(tmp)
#}

ti.ob = get(load(paste(fld,ti.f[k], sep = "/")))
  

Date1 =
data.frame(
samp1= samps$sampleId,
Mo1= samps$min_month,
y1= samps$year
)
Date2 =
data.frame(
samp2= samps$sampleId,
Mo2= samps$min_month,
y2= samps$year
)

ti.ob %<>%
left_join(Date1) %>%
left_join(Date2) %>%
mutate(ydelta = abs(y1-y2))
###
##timeplots
#ti.ob %>%
#  ggplot(aes(
#    x=as.factor(ydelta),
#    y=log(FST),
#    color=pop1
#  )) +
#  geom_boxplot() +
#  facet_wrap(~pop1)->
#  ydelta.plots
#ggsave(ydelta.plots, file = "ydelta.plots.pdf")
#####
#anova(lm(FST ~ pop1*day_diff, data = ti.ob))

####
#### NASA POWER PARAMS
#query_parameters(community = "ag",
#                 temporal_api = "hourly")

###
##ti.ob %>%
##  group_by(pop1) %>%
##  mutate(ymin = min(ydelta)) %>%
##filter(ymin == 0) %>% 
##  filter(ydelta >= 5) %>% 
##  .$pop1 %>% unique -> high.dens.pops
  # = c("Gimenells","Odesa","Munich")
  
####
ti.ob.1y = ti.ob

###########
#### launch
fst.winter.between =
foreach(i = 1:dim(ti.ob.1y)[1],
.combine = "rbind",
.errorhandling = "remove")%do%{

message(i)

tmp = ti.ob.1y[i,]

####
pop.anchor = tmp$pop1
samp1 = tmp$samp1
samp2 = tmp$samp2
date.s1 = as.Date(str_split(samp1, "_")[[1]][5], format = "%Y-%m-%d")
date.s2 = as.Date(str_split(samp2, "_")[[1]][5], format = "%Y-%m-%d")
  
####
nasa.power.sub = nasa.power %>% filter(pop1 == pop.anchor)
nasa.power.sub[which(nasa.power.sub$YYYYMMDD >= min(date.s1,date.s2)
                     &
                       nasa.power.sub$YYYYMMDD <= max(date.s1,date.s2) 
                     ),] %>%
  summarise(T.mean = mean(T2M),
            T.var = var(T2M),
            T.min = min(T2M),
            T.max = max(T2M),
            Tn.below5 = sum(T2M < 5),
            Tn.above32 = sum(T2M > 32),
            ) -> wea.summ

message(wea.summ, sep = "|")

return(data.frame(tmp, wea.summ))

}

###
save(fst.winter.between,
     file = paste("./weather_fst_1y/", "job.weather.fst.", k, ".Rdata", sep ="" )
)
####
