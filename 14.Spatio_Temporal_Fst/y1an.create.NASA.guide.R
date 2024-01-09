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

args = commandArgs(trailingOnly=TRUE)
#k= as.numeric(args[1]) 
##1-923

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

##samps %>%
##  group_by(city) %>%
##  summarize(lat = mean(lat), 
##            long = mean(long)) ->
##  lat.long.samps

###
fld.1y = "/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/outfiles_1y/"

time1y.files = paste(fld.1y,"*",".y1comp",".Rdata", sep = "")
ti1y.f = system(paste("ls ", time1y.files, sep = "" ), intern = T)

#####

ti.1y.ob =
  foreach(i=ti1y.f,
          .combine = "rbind"
  )%do%{
    tmp = get(load(i))
    return(tmp)
  }

#ti.1y.ob = get(load(ti1y.f[k]))

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

ti.1y.ob %<>%
  left_join(Date1) %>%
  left_join(Date2) %>%
  mutate(ydelta = abs(y1-y2))

unique(
unique(ti.1y.ob$samp1),
unique(ti.1y.ob$samp2)
) -> samps.compared

data.frame(
  sampleId= samps.compared
) %>%
  separate(sampleId, remove = F,
           into = c("country",
                    "region",
                    "city",
                    "rep",
                    "date"
                    ), sep = "_") %>%
  filter(!is.na(date)) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) ->
  samps.compared

###########
samps %>%
  group_by(locality) %>%
  summarize(pop1 = unique(city),
            lat = mean(lat),
            long = mean(long)
              ) ->
  cities

samps.compared %>%
  group_by(country,
           region,
           city) %>%
  summarize(minD = min(date),
            maxD = max(date),
            ) %>%
  mutate(locality = 
           paste(country,region,city,
                 sep = "_")) %>%
  left_join(cities) ->
  samps.compared.sums

###############
#### NASA POWER PARAMS
#query_parameters(community = "ag",
#                 temporal_api = "hourly")

#####
Nasapower.winter.range.DEST =
  foreach(i = 1:dim(samps.compared.sums)[1],
          .combine = "rbind")%do%{
            
            message(i)
            
            tmp = samps.compared.sums[i,]
            
            wea.dat <- get_power(
              community = "ag",
              lonlat = c(tmp$long, tmp$lat),
              pars = c("T2M"),
              dates = c(tmp$minD, tmp$maxD),
              temporal_api = "daily",
              time_standard="UTC") 
            
            return(data.frame(tmp, wea.dat))
            
          }

###
save(Nasapower.winter.range.DEST,
     file = "Nasapower.winter.range.DEST.Rdata"
)
####

#wea.dat %>%
#  summarise(T.mean = mean(T2M),
#            T.min = min(T2M),
#            T.max = max(T2M),
#            Tn.below5 = sum(T2M < 5),
#            Tn.above32 = sum(T2M > 32),
#  ) -> wea.summ
#
#message(wea.summ, sep = "|")
