
library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(zoo)
library(adegenet)
library(reshape2)
library("geosphere")
library(gmodels)
library(ggforce)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(forcats)
library(FactoMineR)


load("DEST2.0.model.predictions.Rdata")

out = gims.out.v2
#### plot lat and long
out %>%
  ggplot(
    aes(
      x=`pred.mean$pred.lat`,
      y=real.lat,
      color = predicted.post)
  ) +
  geom_abline(slope = 1) +
  geom_point() ->
  lat.pred

ggsave(lat.pred, file = "lat.pred.pdf")

out %>%
  ggplot(
    aes(
      x=`pred.mean$pred.long`,
      y=real.long,
      color = predicted.post)
  ) +
  geom_abline(slope = 1) +
  geom_point() ->
  long.pred

ggsave(long.pred, file = "long.pred.pdf")

#### Calculate the haversine distance 


  
out.hvd %>%
  ggplot(
    aes(
      x=cont,
      y=(hav_d)
    )
  ) +
  geom_hline(yintercept = 400) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lr")  +
  geom_boxplot() ->
  hav_dist

ggsave(hav_dist, file = "hav_dist.pdf")

out.hvd %>%
  filter(hav_d > 500) %>%
  ggplot(
    aes(
      x=fct_reorder(sampleid, hav_d),
      y=hav_d,
      label=paste(province,predicted.pop,sep = "-->"),
    )
  ) + 
  geom_hline(yintercept = 1000) +
  geom_point() +
  ylim(0, 20000) +
  geom_text(nudge_y = 3000,size = 2) +
  coord_flip() ->
  troubled.samps
  
ggsave(troubled.samps, file = "troubled.samps.pdf")
ggsave(troubled.samps, file = "troubled.samps.png")

