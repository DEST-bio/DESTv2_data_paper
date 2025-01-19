load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE4_Clusters/data_for_reproduction/Figure4_k4.RData")

library(tidyverse)
library(cowplot)
library(pheatmap)

figure1A <- ggplot(a1[a1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1$analysis))
figure1A  

figure1B <- ggplot(b1[b1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
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
figure1B 


figure1_2 <- ggplot(a1_FST_clusters[a1_FST_clusters$type=="Including heterochromatin",], aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
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
figure1_2

heatmap_b2_aut_fgt <- pheatmap(matrix_b2,
                               display_numbers = TRUE,
                               number_color = "black", 
                               fontsize_number = 8,
                               annotation_row = combined_df,
                               annotation_col = combined_df,
                               annotation_colors = ann_colors,
                               main = "Pairwise Fgt - Autosomes")
heatmap_b2_aut_fgt
figure1D <- heatmap_b2_aut_fgt


top<-plot_grid(figure1A+ theme(legend.position="none"),figure1B+ theme(legend.position="none"),align="hv",labels="AUTO") 
bottom<-plot_grid(figure1_2+ theme(legend.position="bottom"),figure1D$gtable,nrow = 2, labels=c("C","D"),rel_heights = c(0.65,1.3))
plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1))
ggsave("fst_global_2.pdf",plot_grid(top,bottom,ncol=1,rel_heights = c(0.25,1)), width = 15, height = 20)
