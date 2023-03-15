###### plot PCA
###### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

###

samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))

####

load("pca.object.Rdata")

pca.object$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) %>%
  filter(!is.na(continent))->
  pca.meta

pca.meta %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = continent
  )) +
  geom_point() +
  theme_bw() ->
  PCA12

ggsave(PCA12, file = "PCA12.pdf", w = 5, h =4)

pca.meta %>%
  filter(continent != "Africa") %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = continent
  )) +
  geom_point() +
  theme_bw() ->
  PCA12na

ggsave(PCA12na, file = "PCA12.noAf.pdf", w = 5, h =4)

##### Run lm with various stats:
##### 
load("/scratch/yey2sn/DEST2_analysis/describe_data/DEST2.0.stats.summary.Rdata")

pca.meta %>%
  left_join(DEST2.0.stats.summary) ->
  pca.meta.stats
  
####
####
cor.test(~Dim.1+lat, data = pca.meta.stats)
cor.test(~Dim.1+long, data = pca.meta.stats)

cor.test(~Dim.1+Missing.data.calc, data = pca.meta.stats)
cor.test(~Dim.1+Neff, data = pca.meta.stats)

summary(lm(Dim.1 ~ lat + long, data = pca.meta.stats))
summary(lm(Dim.2 ~ lat + long, data = pca.meta.stats))

###
###
###
###


#Graph maps
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#### The case of Australia
#####
#### AUSTRALIA
filter(hav.dist.obj.top.pred.city, country == "Australia" ) %>%
  left_join(samps) -> aus.data

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-125.15, 149.99), ylim = c(-48.00, 69.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = pca.meta.stats, 
             aes(x=long, y = lat, fill =Dim.1 ), size = 2.2, shape = 21) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 250 )-> PCA.map
ggsave(PCA.map, file = "PCA.map.pdf", h = 4, w = 9)

