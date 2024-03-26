library(ggplot2)
library(cowplot)
load("figureS10.RData")
figureS1 <- ggplot(a1[a1$type!="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_bw(base_size=16) +
  facet_grid(.~type) +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(a1$analysis))  +
  theme(text=element_text(size=16,color="black"),axis.text = element_text(color="black"))
figureS1


figureS1B <- ggplot(b1[b1$type!="Including\nheterochromatin",], aes(color = chrom, y = mean, x = stat)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se,
                    width = 0.15)) +
  geom_point() + 
  facet_grid(.~type) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_cowplot() +
  ggtitle(unique(b1$analysis)) +
  labs(x = "Statistic", color = "Chromosome", y = "Mean") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#4dbfbc","#ffa500")) +
  ggtitle(unique(b1$analysis))  +
  theme(text=element_text(size=16,color="black"),axis.text = element_text(color="black"))
figureS1B 

ggpubr::ggarrange(figureS1,figureS1B,nrow=2,ncol=1,common.legend = T, legend = "bottom",labels = "AUTO")

ggsave("Figure10.pdf", height = 8, width = 8)
