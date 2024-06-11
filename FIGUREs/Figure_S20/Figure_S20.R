load("Figure_S20.RData")
library(pheatmap)
library(ggplot2)
library(cowplot)
heatmap_b2_aut_fgt <- pheatmap(matrix_b2,
                               display_numbers = TRUE,
                               number_color = "black", 
                               fontsize_number = 8,
                               #annotation_row = combined_df,
                               #annotation_col = combined_df,
                               legend = F,
                               #main = "Pairwise Fgt - Autosomes",
                               annotation_colors = ann_colors)
heatmap_b2_aut_fgt

figure1_2 <- ggplot(a1_FST_clusters[a1_FST_clusters$type=="Including heterochromatin",], aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y = "Mean Fst") +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  #ggtitle(unique(a1_FST_clusters$analysis)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(y=expression("Mean"~italic(F)[ST])) +
  theme(legend.position = "none")
figure1_2

plot_grid(plot_grid(NULL,figure1_2,NULL,ncol=1,rel_heights = c(0.2,1,0.2)), heatmap_b2_aut_fgt$gtable, ncol=2, rel_widths = c(0.7,1), labels="AUTO")

ggsave("Figure_S20.pdf", width = 10, height = 7)
