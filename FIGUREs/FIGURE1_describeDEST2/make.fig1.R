### MAke Figure 1.yes
###

library(data.table)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpattern) ##--> remotes::install_github("coolbutuseless/ggpattern")
library(magick)

setwd("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE1/")
#### samps --->
samps <- fread("code/dest_v2.samps_25Feb2023.csv")

####

#### Contaminantion --->
contam<-get(load("code/Mean.contamination.Final.Rdata"))
#### DuplicateRates --->
duprate<-get(load("code/DuplicateRates.all.Rdata"))
setDT(duprate)
duprate %>%
  dplyr::select(sampleId , pcrdup) ->
  duprate
#### PnPs ---->
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/misc/PnPs/results/classify_pops.txt?token=GHSAT0AAAAAAB7U52CNPDJUDB7GWR4N27GQZCAF6WA -O pnps.predict.txt")
pnps <-fread("code/pnps.predict.txt")
pnps %>%
  filter(Chrom == "genomewide") %>%
  dplyr::select(sampleId = POP, pNpS, private, Status) ->
  pnps.sub
######
pre.dat <-get(load("code/Miss.Cov.PCRdup.sim.joint.Rdata"))
pre.dat %>%
  filter(Var == "Cov") %>%
  dplyr::select(sampleId, Cov = Value) ->
  cov.info

pre.dat %>%
  filter(Var == "Miss") %>%
  dplyr::select(sampleId, Miss = Value) ->
  Miss.info

full_join(contam, duprate) %>%
  full_join(pnps.sub) %>%
  full_join(cov.info) %>%
  full_join(Miss.info) %>%
  full_join(samps) ->
  metadata.seq

####
####

metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) ->metadata.seq

### Panel A
###

samps %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC","dgn") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) %>%
  group_by(set,super.set, continent ) %>%
  summarise(N = n()) %>%
  filter(!is.na(continent)) %>%
  ggplot(aes(set, N, fill = continent)) +
 geom_bar_pattern(stat = "identity",
                 pattern_color = "white",
                 pattern_fill = "black",
                 linewidth = 0.2,
                 aes(pattern = continent))+
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~super.set, ncol = 1,  scales = "free")->
  PanelA

### Panel B
###


samps %>%
  group_by(year, city) %>%
  filter(set != "dest_plus") %>%
  summarize(N = n()) %>%
  .[complete.cases(.),] %>%
  filter(year > 2006) %>%
  filter(N >= 2) %>%
  ggplot(aes(
    y=city,
    x=year,
    fill = N
  )) + geom_point(shape = 22, size = 3.0) +
  scale_fill_gradient2(midpoint = 7.5, mid = "blue", high = "red", low = "grey") +
  theme_classic() + theme(axis.text = element_text(size = 6)) -> PanelB

### Panel C
###

world <- map_data("world")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  ) + theme_classic() -> base_world

samps %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
                               set %in% "dgn" ~ "DGN"
  )) ->samps.map

base_world +
  geom_point(
    data = samps.map,
    aes(x=long,
        y=lat,
        fill = super.set), size = 1.5, shape = 21
  ) +
  scale_fill_brewer(palette = "Set1") -> PanelC

#ggsave(DEST_world, file = "DEST_world.pdf", w = 7, h = 3.5)

### Panel D
metadata.seq %>%
  filter(year >= 2009) %>%
  ggplot(aes(
    y= lat,
    x= jday,
    group=city,
    color=super.set
  )) +
  geom_line() +
  geom_point(size = 1.2) +
  ylim(29,65) +
  xlim(100,350) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
  facet_grid(~year) -> PanelD


#### save panel  PanelA, PanelB, PanelC, PanelD
ggsave(PanelA, file = "PanelA.pdf",
       w= 5, h = 4)

ggsave(PanelB, file = "PanelB.pdf",
       w= 5, h = 7)

ggsave(PanelC, file = "PanelC.pdf",
       w= 6, h = 3)

ggsave(PanelD, file = "PanelD.pdf",
       w= 8, h = 3)
