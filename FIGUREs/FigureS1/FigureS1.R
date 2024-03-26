library(ggplot2)

load("DATA/figS1_df_a.RData")
load("DATA/figS1_df_b.RData")
load("DATA/figS1_df_c.RData")

panelA<- ggplot(figS1_df_a,aes(set, N, fill = continent)) +
  geom_bar(stat = "identity",
                   aes(fill = continent))+
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(legend.position = "top") +
  facet_wrap(~super.set, ncol = 1,  scales = "free") +
  labs(y="",fill="Continent",x="")  +
  theme(text=element_text(size=16,color="black"),axis.text = element_text(color="black")) +
  scale_fill_brewer(palette = "Spectral") 
  
panelA

panelB <-  ggplot(figS1_df_b,aes(
    y=city,
    x=year,
    fill = N
  )) + geom_point(shape = 22, size = 3.0) +
    scale_fill_gradient2(midpoint = 7.5, mid = "#3288bd", high = "#d53e4f", low = "lightgrey") +
    theme_bw(base_size = 16)  +
  theme(text=element_text(size=16,color="black"),axis.text = element_text(color="black"), legend.position = "right") +  labs(x="",y="") 
panelB

panelC <- 
  ggplot() + 
  geom_point(
    data=figS1_df_c,
    aes(
      x=rank.i,
      y=(value),
      #color = set,
      #shape = set
    ), size = 2, alpha = 0.35) +
  #geom_point(data=plot.data[collapse == "Yes"],
  #           aes(
  #             x=rank.i,
  #             y=(value),
  #             #shape = collapse,
  #             #fill = set,
  #           ), shape = 21, size = 2, alpha = 0.3, color = "black") +
  facet_grid(variable~set, scales = "free") +
  theme_bw(base_size = 16) +
  ylab("") +
  xlab("Sample") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=16,color="black"),axis.text = element_text(color="black"))
panelC

cowplot::plot_grid(panelA,panelB,panelC,ncol = 1,labels="AUTO",align="l",rel_heights = c(0.6,1.4,1))
ggsave("figureS1.pdf",width = 12, height = 22)
