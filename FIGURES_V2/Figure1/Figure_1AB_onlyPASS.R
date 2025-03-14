# Libraries
library(data.table)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpattern) 
library(magick)
library(metR)
library(plyr)
library(ggthemes)
library(ggpol)
library(scales)

setwd("~/Desktop/DEST/Figure1")

##### A) Map
#### samples file
world <- map_data("world")

last_metadata <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")
last_metadata <- last_metadata %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
                               set %in% "dgn" ~ "DGN"
  ))

a <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    fill = "#cccccc", color = "darkgray", linewidth = 0.1
  ) +
  theme_classic(base_size = 14) +
  scale_x_longitude(breaks = seq(-160,180,40)) +
  scale_y_latitude(breaks = seq(-40,100,20),limits=c(-61,76)) +
  coord_sf() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_fill_manual(values = c("#f27f65","#8c89c1","#a5c9cc")) +
  scale_color_manual(values = c("#D05438","#615E93","#6F9FA3")) +
  labs(fill="Set", color="Set") +
  theme(legend.position = "bottom") + 
  geom_point(
    data = last_metadata[last_metadata$Recommendation=="Pass",],
    aes(x=long,
        y=lat,
        fill = super.set, color = super.set), size = 2.5, shape = 21
  )
a


# only with pass

last_metadata.df <- last_metadata %>%
  filter(set != "dgn") %>%
  filter(year >= 2009)

b <- ggplot(last_metadata.df[last_metadata.df$Recommendation=="Pass",], aes(
  y= lat,
  x= jday,
  group=city,
  color=super.set,fill=super.set
)) +
  geom_line(alpha=0.5, show.legend = F) +
  geom_point(show.legend = F,
             aes(fill = super.set, color = super.set), size = 1.5, shape = 21
  ) +
  scale_y_latitude(name = "Latitude") +
  scale_x_continuous(limits=c(100,350),breaks = seq(100,400,100)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  labs(x="Julian day", y = "Latitude")+
  facet_grid(~year)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_fill_manual(values = c("#f27f65","#8c89c1","#a5c9cc")) +
  scale_color_manual(values = c("#D05438","#615E93","#6F9FA3")) +
  labs(fill="Set", color="Set") +
  theme(legend.position = "bottom")
b

library(cowplot)
ggsave("figure1_v4.pdf",plot=plot_grid(a,b,align="ltrb",labels=c("A","B"),ncol=1,rel_heights = c(1,0.6)), height = 8, width = 8)
