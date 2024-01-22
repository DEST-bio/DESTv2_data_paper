library(ggplot2)
library(data.table)
library(dplyr)
library(DT)
library(eulerr)
library(cowplot)
library(stringr)
library(tidyr)
library(ggrepel)
library(UpSetR)
library(ggthemes)
library(tibble)
library(pheatmap)
library(scales)

load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FSG_plots/spatialFST_figure.RData")

# main figure spatial fst
figureA <- ggplot(a1[a1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1$analysis))
figureA  

figureB <- ggplot(b1[b1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  ggtitle(unique(b1$analysis)) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(b1$analysis))
figureB   


figureC <- ggplot(a1_FST_clusters[a1_FST_clusters$type=="Including heterochromatin",], aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y = "Mean Fst") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1_FST_clusters$analysis)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
figureC

figureD <- pheatmap(matrix_b2,
                               display_numbers = TRUE,
                               number_color = "black", 
                               fontsize_number = 8,
                               annotation_row = combined_df,
                               annotation_col = combined_df,
                               annotation_colors = ann_colors,
                               main = "Pairwise Fgt - Autosomes")
figureD


top<-plot_grid(figureA+ theme(legend.position="none"),figureB+ theme(legend.position="none"),align="hv",labels="AUTO") 
bottom<-plot_grid(figureC+ theme(legend.position="bottom"),figureD$gtable,nrow = 2, labels=c("C","D"),rel_heights = c(0.65,1.3))
plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1))
ggsave("fst_global.pdf",plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1)), width = 15, height = 20)
ggsave("fst_global.png",plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1)), width = 15, height = 20)


figureSA <- ggplot(a1[a1$type!="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  facet_grid(.~type) +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1$analysis))
figureSA


figureSB <- ggplot(b1[b1$type!="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  facet_grid(.~type) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  ggtitle(unique(b1$analysis)) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(b1$analysis))
figureSB 


figureSC <- ggplot(a1_FST_clusters[a1_FST_clusters$type!="Including heterochromatin",], aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y = "Mean Fst") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1_FST_clusters$analysis)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(.~type)
figureSC

top<-plot_grid(figureSA+ theme(legend.position="none"),figureSB+ theme(legend.position="none"),align="hv",labels="AUTO") 
bottom<-plot_grid(figureSC,labels=c("C"))
plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1))

