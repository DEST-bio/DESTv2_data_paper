
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(ggExtra)
library(foreach)
library(data.table)
library(factoextra)
library(cowplot)
library(pheatmap)

#load geoanalysis packages
library(rnaturalearth)
library(rnaturalearthdata)

dataset <- read.table("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE2_Invs/FullInvDestv2.txt", header = T, na.string = "NA")
# dataset=read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/EurInv.txt",header=T,na.string="NA")

dataset <- na.omit(dataset)
attach(dataset)
dataset %>% melt(id = c("Latitude",   "Longitude", "DEST") ) -> dataset.m

# the code below creates a "base world" figure... needed for all maps
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1
  ) + theme_classic() -> base_world

### plot
base_world + 
  coord_sf(xlim =  c(-135, 160), ylim = c(-62, 68.00), expand = FALSE) + 
  geom_point(
    data = dataset.m,
    aes(x=Longitude,
        y=Latitude,
        fill = value,
        shape= as.factor(DEST)
        ), 
    size = 2.3
  ) + scale_shape_manual(values = 21:22) + 
  scale_fill_gradientn(
    name = "Inversion Frequency",
    colours = c("white", "yellow", "red", "purple", "blue")) +
  facet_wrap(~variable, ncol = 2) -> world_plot

ggsave("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE2_Invs/Fig2.new.pdf",
    world_plot,
    width = 12,
    height = 6
)
