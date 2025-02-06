###### plot PCA
###### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)

samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

###
### GET DATA
#data from "https://github.com/DEST-bio/DESTv2_data_paper/raw/refs/heads/main/8.PCA_PopStruct/data/PCA.results.df.Rdata"
all.dat.pca <- get(load("/Users/jcnunez/Downloads/PCA.results.df.Rdata"))
setDT(all.dat.pca)

### Plot ->
all.dat.pca %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color=continent
  )) +
  geom_point(size = 0.8) +
  theme_classic() +
  scale_color_brewer(palette = "Spectral" ) +
  theme(legend.position = "top") +
  facet_grid(~case) -> PC12.chrs
ggsave(PC12.chrs, file = "PC12.chrs.pdf", w = 8, h = 3.0)

all.dat.pca %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.3,
    color=continent
  )) +
  geom_point(size = 0.8) +
  theme_classic() +
  scale_color_brewer(palette = "Spectral" ) +
  theme(legend.position = "top") +
  facet_grid(~case) -> PC13.chrs
ggsave(PC13.chrs, file = "PC13.chrs.pdf", w = 8, h = 3.0)

#### see perc explained
#### 
perc <- get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE2_PCA/data_for_reproduction/pca.var.exp.Rdata"))
perc = do.call(rbind, perc)
setDT(perc)
perc %>%
  filter(PC %in% c("comp 1","comp 2","comp 3")) %>% 
  filter(PC %in% c("comp 3"))

### Some correlations
### 
all.dat.pca %>%
  filter(case == "all") %>%
  filter(continent != "Africa") %>%
  filter(country != "China") %>%
  dplyr::select(Dim.1,Dim.2,Dim.3,long,continent) %>%
  melt(id = c("long","continent") ) %>% 
  ggplot(aes(
    x=long,
    y=value,
    fill=continent,
  )) +
  geom_point(shape = 21) +
  geom_smooth(method = "lm", color = "grey20", size = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral" ) +
  theme(legend.position = "none") +
  facet_grid(variable~continent, scales = "free_x") -> Pc1.long
ggsave(Pc1.long, file = "Pcs.long.pdf", w = 8, h = 2.7)

all.dat.pca %>%
  filter(case == "all") %>%
  filter(continent != "Africa") %>%
  filter(country != "China") %>%
  dplyr::select(Dim.1,Dim.2,Dim.3,lat,continent) %>%
  melt(id = c("lat","continent") ) %>% 
  ggplot(aes(
    x=lat,
    y=value,
    fill=continent,
  )) +
  geom_point(shape = 21) +
  geom_smooth(method = "lm", color = "grey20", size = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral" ) +
  theme(legend.position = "none") +
  facet_grid(variable~continent, scales = "free_x") -> Pc1.lat
ggsave(Pc1.lat, file = "Pcs.lat.pdf", w = 8, h = 2.7)


####
####
cor.test(~Dim.1+lat, data = filter(all.dat.pca, case == "all", continent == "Europe"))
cor.test(~Dim.2+lat, data = filter(all.dat.pca, case == "all", continent == "Europe"))
cor.test(~Dim.3+lat, data = filter(all.dat.pca, case == "all", continent == "Europe"))

cor.test(~Dim.1+long, data = filter(all.dat.pca, case == "all", continent == "Europe"))
cor.test(~Dim.2+long, data = filter(all.dat.pca, case == "all", continent == "Europe"))
cor.test(~Dim.3+long, data = filter(all.dat.pca, case == "all", continent == "Europe"))

cor.test(~Dim.1+lat, data = filter(all.dat.pca, case == "all", continent == "North_America"))
cor.test(~Dim.2+lat, data = filter(all.dat.pca, case == "all", continent == "North_America"))
cor.test(~Dim.3+lat, data = filter(all.dat.pca, case == "all", continent == "North_America"))

cor.test(~Dim.1+long, data = filter(all.dat.pca, case == "all", continent == "North_America"))
cor.test(~Dim.2+long, data = filter(all.dat.pca, case == "all", continent == "North_America"))
cor.test(~Dim.3+long, data = filter(all.dat.pca, case == "all", continent == "North_America"))

#### Compare PCs 
dim = "Dim.3"
all.dat.pca %>% 
  filter(case %in% c("all","2L")) %>%
  dcast(sampleId~case, value.var = dim) %>%
  cor.test(~`2L`+all, data = .)
all.dat.pca %>% 
  filter(case %in% c("all","2R")) %>%
  dcast(sampleId~case, value.var = dim) %>%
  cor.test(~`2R`+all, data = .)
all.dat.pca %>% 
  filter(case %in% c("all","3L")) %>%
  dcast(sampleId~case, value.var = dim) %>%
  cor.test(~`3L`+all, data = .)
all.dat.pca %>% 
  filter(case %in% c("all","3R")) %>%
  dcast(sampleId~case, value.var = dim) %>%
  cor.test(~`3R`+all, data = .)

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


##### MBE REVISION -- PCA and Metadata
#####
##### MBE REVISION -- PCA and Metadata
#####
##### MBE REVISION -- PCA and Metadata
#####
##### MBE REVISION -- PCA and Metadata
#####
##### MBE REVISION -- PCA and Metadata
#####

samps %>%
  #filter(city %in% c("Yesiloz","Odesa","Charlottesville")) %>%
  mutate(sampling_curated=
  case_when(
    sampling_strategy %in% c("Sweep netting", "Nets", "Sweep Netting") ~ "Netting",
    sampling_strategy %in% c("Aspirating", "Aspirator") ~ "Aspirator",
    sampling_strategy %in% c("Fly trap") ~ "Fly_trap",
    TRUE ~ "multiple"
  )) %>%
  filter(fruit_type_curated %in% c("Grape","Apple","Watermelon","Peach","Pear","Cherry","Beer", "Banana")) ->
  samps.cur

all.dat.pca %>%
  filter(case == "all") %>%
  left_join(samps.cur) %>%
  #filter(city %in% c("Yesiloz","Odesa","Charlottesville")) %>%
  mutate(Neff = (Cov*nFlies)/(Cov+nFlies-1) ) %>%
  filter(!is.na(Neff)) ->
  pca.all




ANOVAS =
foreach(i=2:4, .combine = "rbind")%do%{
  
  lm(Dim.1 ~ Neff + lat*long + as.factor(year) + jday*fruit_type_curated*sampling_curated, 
     data = filter(pca.all, cluster2.0_k4 == i)) %>% anova -> an1
  lm(Dim.2 ~ Neff + lat*long + as.factor(year) + jday*fruit_type_curated*sampling_curated,
     data = filter(pca.all, cluster2.0_k4 == i)) %>% anova -> an2
  lm(Dim.3 ~ Neff + lat*long + as.factor(year) + jday*fruit_type_curated*sampling_curated, 
     data = filter(pca.all, cluster2.0_k4 == i)) %>% anova -> an3
  
rbind(
data.frame(cluster = i, PC = 1, an1), 
data.frame(cluster = i, PC = 2, an2), 
data.frame(cluster = i, PC = 3, an3) 
)
}

write.table(ANOVAS, file = "ANOVAS.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


pca.all %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=(fruit_type_curated),
    shape = (sampling_curated)
  )) + geom_point(size = 3)  + theme_bw() +
  scale_shape_manual(values = 21:24) +
  facet_wrap(~cluster2.0_k4, scales = "free", nrow = 1) -> fruits
ggsave(fruits, file = "fruits.pdf", w = 9, h = 3)

pca.all %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=jday,
    shape = as.factor(sampling_curated)
  )) + geom_point(size = 3, ) + theme_bw() +
  scale_shape_manual(values = 21:24) +
  scale_fill_gradient2(midpoint = 250, high = "blue", low = "red") +
  facet_wrap(~cluster2.0_k4, scales = "free", nrow = 1) -> traps_time

ggsave(traps_time, file = "traps_time.pdf", w = 9, h = 3)
