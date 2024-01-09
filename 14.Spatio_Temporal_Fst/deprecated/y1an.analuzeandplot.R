##### 1-y analyses


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

##
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

