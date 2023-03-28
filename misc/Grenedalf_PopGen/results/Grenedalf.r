

library(tidyverse)
library(cowplot)

DATA=read.table("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf_poolsnp.summary",header=T,sep="\t")

DATA.gw.theta_pi.p<-na.omit(DATA) %>%
    filter(DATA$Chrom=="GenomeWide"& DATA$Stat!="snp_count")

P.pool <- ggplot(DATA.gw.theta_pi.p,aes(x=Continent,y=as.numeric(Value)))+
    geom_boxplot()+
    ggtitle("PoolSNP")+
    facet_wrap(.~Stat,scales="free_y")+
    ylab("PopGen Statistic")+
    theme_bw()

DATA=read.table("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf_snape.summary",header=T,sep="\t")

DATA.gw.theta_pi.s<-na.omit(DATA) %>%
    filter(DATA$Chrom=="GenomeWide"& DATA$Stat!="snp_count")

P.snape <- ggplot(DATA.gw.theta_pi.s,aes(x=Continent,y=as.numeric(Value)))+
    geom_boxplot()+
    ggtitle("SNAPE")+
    facet_wrap(.~Stat,scales="free_y")+
    ylab("PopGen Statistic")+
    theme_bw()

P<-plot_grid(P.pool,P.snape,
    nrow=2)

DATA.gw.theta_pi.s$Caller <- rep("SNAPE",nrow(DATA.gw.theta_pi.s))
DATA.gw.theta_pi.p$Caller <- rep("PoolSNP",nrow(DATA.gw.theta_pi.p))
DATA.new<-spread(rbind(DATA.gw.theta_pi.s,DATA.gw.theta_pi.p),Caller,Value)


P <- ggplot(DATA.new,aes(x=SNAPE,y=PoolSNP,col=Stat))+
    geom_point()+
     facet_wrap(.~Stat,scales="free")+
    geom_smooth(
    method = "lm",
    formula = y ~ x)+
    theme_bw()


ggsave("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf_corr.pdf",
    P,
    width=16,
    height=8)

ggsave("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf_corr.png",
    P,
    width=16,
    height=8)



