#### Final figure

#### Collect and plot admixture proportions 
library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(poolfstat)
library(FactoMineR)
require(foreach)
require(gmodels)
library(rnaturalearth)
library(rnaturalearthdata)
require(DescTools)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


datin <- get(load("/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/ancestry_Estimates_linear_moments.Rdata")) 

datin %>%
  ggplot(aes(
    x=lat,
    y=1-ancestry,
    shape = continent
  )) + 
  geom_smooth(method = "lm", se = T) +
  geom_point(size = 2.1, fill = "grey") +
  scale_shape_manual(values = 21:23) +
  facet_grid(~continent, scales = "free_x") +
  theme_bw() -> momemts_ests

ggsave(momemts_ests, file = "momemts_ests.pdf",
       w= 8, h = 2.5)

####
summary(lm(c(1-ancestry) ~ lat, data = filter(datin, continent == "North_America") ))
summary(lm(c(1-ancestry) ~ lat, data = filter(datin, continent == "Oceania") ))
summary(lm(c(1-ancestry) ~ lat, data = filter(datin, continent == "South_America") ))

summary(lm(c(1-ancestry) ~ long, data = filter(datin, continent == "North_America") ))


########
datin %>%
  filter(continent == "North_America") %>%
  group_by(lat,long) %>%
  summarize(m.est = mean(1-ancestry)) -> NAMAncestry

## world plot
ggplot(data = world) +
  geom_sf(fill= "grey70") +
  coord_sf(xlim = c(-50, -135), ylim = c(15,50), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = NAMAncestry , 
              aes(
                x=long,
                y=lat,
                fill = m.est,
              ), alpha = 0.9, size = 3.5, shape = 21, color = "white") + 
  scale_fill_viridis_c()  ->
  WordAFR.admix

ggsave( WordAFR.admix, file = "moments.WordAFR.admix.pdf", w = 9, h  =3)

####
datin %>%
  ggplot(aes(
    x=1-ancestry,
    y=linearEST,
    shape = continent
  )) +  geom_abline(slope = 1) +
  geom_point(fill = "grey", size = 2.1) +
  scale_shape_manual(values = 21:23) +
  ylim(0,0.4) + xlim(0,0.4) +
  theme_bw() -> corr_methods

ggsave(corr_methods, file = "corr_methods.pdf",
       w= 4.5, h = 3)

###
cor.test(~c(1-ancestry)+linearEST, data = datin)

### pattern
datin %>%
  filter(continent == "North_America") %>%
  mutate(est_diff = c(1-ancestry)+linearEST) %>%
  ggplot(aes(
    x=c(1-ancestry),
    y=est_diff,
  )) + 
  geom_abline(slope = 1) +
  geom_smooth(method = "lm", se = T) +
  geom_point(shape = 21, size = 2.1, fill = "grey") +
  facet_grid(~continent, scales = "free_x") +
  xlim(0.05,0.5) + ylim(0.05,0.5)+
  theme_bw() -> diff_plot.MOLI

ggsave(diff_plot.MOLI, file = "diff_plot.MOLI.pdf",
       w= 4, h = 4)

cc1 <- CCC(c(1-(datin$ancestry)), datin$linearEST, na.rm = T)

###
datin %>%
  group_by(cluster2.0_k8) %>%
  ggplot(aes(
    x=as.character(cluster2.0_k8),
    y=1-ancestry
  )) + geom_boxplot() ->
  box_ance

ggsave(box_ance, file = "box_ance.pdf")

