### Leverage analysis


### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
library(vroom)
library(reshape2)

meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv"
samps <- fread(meta_git)
setDT(samps)

samps %>%
  dplyr::select(focal=sampleId, focal_contnent = continent) ->
  samps.focal

load("linear.admix.dat.Rdata")
linear.admix.dat.filters %>%
  filter(filter == "noINV") %>%
  filter(source_pop == "AFRICA") %>%
  filter(admix.set == "N.America") ->
  nAm.data

model <- lm(mean.est ~ lat , data=nAm.data)
#view model summary
summary(model)
hats <- as.data.frame(hatvalues(model))

#display leverage stats for each observation
nAm.data %<>%
  as.data.frame() %>%
  mutate(leverage=hats$hatvalues)

nAm.data %>%
  ggplot(aes(
    x=lat,
    y=leverage
  )) + geom_point() + geom_hline(yintercept = 2) ->
  leverage.analisis

ggsave(leverage.analisis, file = "leverage.analisis.png",
       h=4, w=4)
###
model2 <- lm(mean.est ~ lat , data=filter(nAm.data, lat > 20))
#view model summary
summary(model2)


