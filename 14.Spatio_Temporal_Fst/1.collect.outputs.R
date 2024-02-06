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

### wheather
#nasapower <- get(load("/netfiles/nunezlab/Drosophila_resources/NASA_power_weather/DEST2.0/nasaPower.allpops.Rdata"))

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

samps %>%
group_by(city) %>%
summarize(lat = mean(lat), 
long = mean(long)) ->
lat.long.samps

###
fld = "/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/outfiles/"

#space.files = paste(fld,"*",".space",".Rdata", sep = "")
#sp.f = system(paste("ls ", space.files, sep = "" ), intern = T)

####
#sp.ob =
#foreach(i=sp.f,
#.combine = "rbind"
#)%do%{
#tmp = get(load(i))
#return(tmp)
#}

###
time.files = paste(fld,"*",".time",".Rdata", sep = "")
ti.f = system(paste("ls ", time.files, sep = "" ), intern = T)

####
ti.ob =
foreach(i=ti.f,
.combine = "rbind"
)%do%{
tmp = get(load(i))
return(tmp)
}

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
filter(pop1 != "Vully") %>%
left_join(Date1) %>%
left_join(Date2) %>%
mutate(ydelta = abs(y1-y2))
###
##timeplots
ti.ob %>%
  ggplot(aes(
    x=as.factor(ydelta),
    y=log(FST),
    color=pop1
  )) +
  geom_boxplot() +
  facet_wrap(~pop1)->
  ydelta.plots
ggsave(ydelta.plots, file = "ydelta.plots.pdf")
####
anova(lm(FST ~ pop1*day_diff, data = ti.ob))

####
#### NASA POWER PARAMS
query_parameters(community = "ag",
                 temporal_api = "hourly")

###
##ti.ob %>%
##  group_by(pop1) %>%
##  mutate(ymin = min(ydelta)) %>%
##filter(ymin == 0) %>% 
##  filter(ydelta >= 5) %>% 
##  .$pop1 %>% unique -> high.dens.pops
  # = c("Gimenells","Odesa","Munich")
  
####
ti.ob %>%
filter(pop1 != "Vully") ->
ti.ob.1y

#### launch
fst.winter.between =
foreach(i = 1:dim(ti.ob.1y)[1],
.combine = "rbind")%do%{

message(i)

tmp = ti.ob.1y[i,]

loc = tmp$samp1
site= tmp$pop1

lat.loc=samps$lat[which(samps$sampleId == loc)]
long.loc=samps$long[which(samps$sampleId == loc)]

target_year1 = min(tmp$y1,tmp$y2)
target_year2 = max(tmp$y1,tmp$y2)

if(target_year1 < target_year2){
  trgMo1 = tmp$Mo1
  trgMo2 = tmp$Mo2
  
  trgDa1= str_split(str_split(tmp$samp1, "_")[[1]][5], "-")[[1]][3]
  trgDa2= str_split(str_split(tmp$samp2, "_")[[1]][5], "-")[[1]][3]
    
} else{
  trgMo1 = tmp$Mo2
  trgMo2 = tmp$Mo1
  trgDa1= str_split(str_split(tmp$samp2, "_")[[1]][5], "-")[[1]][3]
  trgDa2= str_split(str_split(tmp$samp1, "_")[[1]][5], "-")[[1]][3]
}

if(nchar(trgMo1) == 1){
  trgMo1 = paste("0",trgMo1, sep ="")
} else if(nchar(trgMo1) > 1){trgMo1=trgMo1}

if(nchar(trgMo2) == 1){
  trgMo2 = paste("0",trgMo2, sep ="")
} else if(nchar(trgMo2) > 1){trgMo2=trgMo2}


wea.dat <- get_power(
  community = "ag",
  lonlat = c(long.loc, lat.loc),
  pars = c("T2M"),
  dates = c(paste(target_year1,trgMo1 ,trgDa1, sep="-"), 
            paste(target_year2,trgMo2 ,trgDa2, sep="-")),
  temporal_api = "daily",
  time_standard="UTC") 

wea.dat %>%
  summarise(T.mean = mean(T2M),
            T.min = min(T2M),
            T.max = max(T2M),
            Tn.below5 = sum(T2M < 5),
            Tn.above32 = sum(T2M > 32),
            ) -> wea.summ

message(wea.summ, sep = "|")

return(data.frame(tmp, wea.summ))

}

save(fst.winter.between,
     file = "fst.winter.newTimePops.Rdata"
)
