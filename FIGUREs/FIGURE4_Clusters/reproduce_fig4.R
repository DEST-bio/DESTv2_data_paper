#### Reconstruct Figure 4 as. of Feb 26, 2024
#### Jcbn

#load quality of life
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

####
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

#### Panel A/D -- mapf of K=4/8

# the code below creates a "base world" figure... needed for all maps
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1
  ) + theme_classic() -> base_world


### the code below extracts the PCs from the PCA object
base_world + 
  geom_point(
    data = samps[complete.cases(samps$cluster2.0_k8),],
    aes(x=long,
        y=lat,
        fill = as.factor(cluster2.0_k4)), size = 2.3, shape = 21
  ) +
  scale_fill_brewer(palette = "Set1") -> cluster.plot.k4
ggsave(cluster.plot.k4, file = "cluster.plot.k4.pdf", w = 7, h = 3.5)

base_world + 
  geom_point(
    data = samps[complete.cases(samps$cluster2.0_k8),],
    aes(x=long,
        y=lat,
        fill = as.factor(cluster2.0_k8)), size = 2.3, shape = 21
  ) +
  scale_fill_brewer(palette = "Set2") -> cluster.plot.k8
ggsave(cluster.plot.k8, file = "cluster.plot.k8.pdf", w = 7, h = 3.5)

### Create B and C
### these are the zoom in plots
#### Zoom on Europe
world2 <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world2) +
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-12, 41.00), ylim = c(32.00, 63.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(
             data = samps[complete.cases(samps$cluster2.0_k8),],
             aes(x=long,
                 y=lat,
                 fill = as.factor(cluster2.0_k8)), size = 2.3, shape = 21) +
  xlab("Lon") + ylab("Lat") + 
  ggtitle("Europe under K=8") + 
  theme(legend.position = "none") + 
  scale_shape_manual(values = c(21,22,23,24)) -> Suture_ZoneEU
ggsave(Suture_ZoneEU, file = "Suture_ZoneEU.pdf",  w = 7, h = 3.5)

### Plot North America
ggplot(data = world2) +
  geom_sf(fill= "lightgray") +
  coord_sf(xlim =  c(-125.15, -55), ylim = c(-10, 50.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(
    data = samps[complete.cases(samps$cluster2.0_k8),],
    aes(x=long,
        y=lat,
        fill = as.factor(cluster2.0_k8)), size = 2.3, shape = 21) +
  xlab("Lon") + ylab("Lat") + 
  ggtitle("North America under K=8") + 
  theme(legend.position = "none") + 
  scale_shape_manual(values = c(21,22,23,24)) -> Suture_ZoneAM
ggsave(Suture_ZoneAM, file = "Suture_ZoneAM.pdf",  w = 7, h = 3.5)


#### Make panels E-G --> Linear admixture plots
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE4_Clusters/data_for_reproduction/linear.admix.dat.Guinea.EUW.Rdata")

linear.admix.dat.filters %>%
  filter(source_pop == "AFRICA") %>%
  filter(filter %in% c("noINV","All")) %>%
  ggplot(aes(
    x=lat,
    y=mean.est,
    ymin = mean.est - sd.est,
    ymax = mean.est + sd.est,
    color = source_pop
  )) + 
  geom_smooth(method = "lm",aes(linetype = filter)) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2.1, 
             fill = "grey",
             color = "black", aes(shape = filter)) +
  theme_bw() +
  scale_shape_manual(values = 21:22) +
  facet_grid(.~admix.set, scales = "free_x")->
  plot.admix.flt

ggsave(plot.admix.flt, file = "plot.admix.flt.pdf", w = 10, h  =3)

### Panels H-J
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE4_Clusters/data_for_reproduction/spatialFST_figure.RData")

figureH <- ggplot(a1[a1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point(size = 3) + 
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1$analysis))
figureH  
ggsave(figureH, file = "figureH.pdf", h = 2.5, w = 3)

figureI <- ggplot(b1[b1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point(size = 3) + 
  #facet_grid(.~type) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  ggtitle(unique(b1$analysis)) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(b1$analysis))
figureI   
ggsave(figureI, file = "figureI.pdf", h = 2.5, w = 3.5)


figureK <- ggplot(a1_FST_clusters[a1_FST_clusters$type=="Including heterochromatin",], aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point(size = 3) + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y = "Mean Fst") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1_FST_clusters$analysis)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
figureK
ggsave(figureK, file = "figureK.pdf", h = 2.5, w = 3.5)

figureJ <- pheatmap(matrix_b2,
                    display_numbers = TRUE,
                    number_color = "black", 
                    fontsize_number = 8,
                    annotation_row = combined_df,
                    annotation_col = combined_df,
                    annotation_colors = ann_colors,
                    main = "Pairwise Fgt - Autosomes")
figureJ
ggsave(figureJ, file = "figureJ.pdf", h = 6, w = 9)
