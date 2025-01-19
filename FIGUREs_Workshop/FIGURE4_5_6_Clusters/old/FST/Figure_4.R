load("Figure_4.RData")
library(pheatmap)
library(ggplot2)
library(cowplot)

matrix_b2_tr<-matrix_b2
matrix_b2_tr[lower.tri(matrix_b2_tr)] <- NA

plotH <- pheatmap(matrix_b2_tr,
                               display_numbers = TRUE,
                               number_color = "black", 
                               fontsize_number = 12, fontsize = 12,cluster_rows = F, cluster_cols = F, na_col="white",border_color ="white", main=expression(italic(F)[GT]))
plotH

plotJ <- ggplot(a1[a1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  labs(x = "", color = "Chromosome", y = expression("Mean"~italic(F)[ST])) +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  theme(legend.position = "none")  +
  scale_x_discrete(labels="DEST 2.0")
  #ggtitle(unique(a1$analysis))
plotJ  




plotK <- ggplot(b1[b1$type=="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  #ggtitle(unique(b1$analysis)) +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c(expression(italic(F)[GT]), expression(italic(F)[SG]), expression(italic(F)[ST])))
  #ggtitle(unique(b1$analysis))
plotK  



plotL <- ggplot(a1_FST_clusters, aes(color = chrom, y = mean, x = reorder(cluster,mean))) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  #facet_grid(.~type) +
  labs(x = "", color = "Chromosome", y =expression("Mean"~italic(F)[ST])) +
  theme_cowplot() +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  #ggtitle(unique(a1_FST_clusters$analysis)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plotL

plot_grid(plotJ,plotK,plotL,ncol = 3, rel_widths = c(0.5,0.5,1), align = "h")
ggsave("4JKL.pdf", height = 4, width = 10)