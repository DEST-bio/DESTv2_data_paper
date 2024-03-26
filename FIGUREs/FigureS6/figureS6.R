library(ggplot2)
library(scales)
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
pca1<-load("PCA.genome.euch.w200k.rel.n.ind.clusters.nolabels.RData")
pca1<-get(pca1[[1]])
pca2<-load("PCA.genome.euch.w200k.rec.n.ind.clusters.nolabels.RData")
pca2<-get(pca2[[1]])

plot_3R<-load("3R.RData")
plot_3R<-get(plot_3R[[1]])
plot_3R_df <- plot_3R$data


plot_3L<-load("3L.RData")
plot_3L<-get(plot_3L[[1]])

plot_2R<-load("2R.RData")
plot_2R<-get(plot_2R[[1]])

plot_2L<-load("2L.RData")
plot_2L<-get(plot_2L[[1]])

plot_X<-load("X.RData")
plot_X<-get(plot_X[[1]])

plot_2L <- plot_2L + labs(y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_segment(aes(x = 2225744, xend = 13154180, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 13154180, xend = 2225744, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 6955744, y = 0.000000035, label = "In(2L)t",fontface = 'italic',fill="white") 


plot_2R <- plot_2R + labs(y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_segment(aes(x = 15391154, xend = 20276334, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 20276334, xend = 15391154, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 17833744, y = 0.000000035, label = "In(2R)NS",fontface = 'italic',fill="white") 


plot_3L <- plot_3L + labs(y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_segment(aes(x = 3173046, xend = 16308841, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 16308841, xend = 3173046, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 9740944, y = 0.000000035, label = "In(3L)P",fontface = 'italic',fill="white") 


plot_3R<-plot_3R + labs(y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_segment(aes(x = 21406917, xend = 29031297, y = 0.000000037, yend = 0.000000037),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 29031297, xend = 21406917, y = 0.000000037, yend = 0.000000037),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 25219107, y = 0.000000037, label = "In(3R)Mo",fontface = 'italic',fill="white")  +
  geom_segment(aes(x = 16432209, xend = 24744010, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 24744010, xend = 16432209, y = 0.000000035, yend = 0.000000035),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 20588110, y = 0.000000035, label = "In(3R)P",fontface = 'italic',fill="white")  +
  geom_segment(aes(x = 20096867, xend = 32079331, y = 0.000000033, yend = 0.000000033),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 32079331, xend = 20096867, y = 0.000000033, yend = 0.000000033),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 26088099, y = 0.000000033, label = "In(3R)C",fontface = 'italic',fill="white") +
  geom_segment(aes(x = 11750567, xend = 26140370, y = 0.000000031, yend = 0.000000031),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 26140370, xend = 11750567, y = 0.000000031, yend = 0.000000031),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 18945468, y = 0.000000031, label = "In(3R)K",fontface = 'italic',fill="white") 


plot_X <- plot_X + labs(y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_segment(aes(x = 13625736, xend = 17828912, y = 0.000000006, yend = 0.000000006),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 17828912, xend = 13625736, y = 0.000000006, yend = 0.000000006),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 15727324, y = 0.000000006, label = "In1A",fontface = 'italic',fill="white")  +
  geom_segment(aes(x = 17828912, xend = 19593711, y = 0.000000005, yend = 0.000000005),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1) +
  geom_segment(aes(x = 19593711, xend = 17828912, y = 0.000000005, yend = 0.000000005),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt", color = "firebrick4", size = 0.1)  +
  annotate("label", x = 18711312, y = 0.000000005, label = "InBe",fontface = 'italic',fill="white") 

panel_A <- cowplot::plot_grid(plot_2L,plot_2R,plot_3L,plot_3R,plot_X,ncol=1)

panel_B <- cowplot::plot_grid(pca1,pca2,NULL,NULL,ncol=1)

cowplot::plot_grid(panel_A,panel_B,ncol=2,labels="AUTO",rel_widths = c(1,0.8))

ggsave("FigureS6.pdf", height = 25, width = 15)
