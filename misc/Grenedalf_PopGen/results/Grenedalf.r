

library(tidyverse)

DATA=read.table("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf.summary",header=T,sep="\t")

DATA.gw.theta_pi<-na.omit(DATA) %>%
    filter(DATA$Chrom=="GenomeWide"& DATA$Stat!="snp_count")

P <- ggplot(DATA.gw.theta_pi,aes(x=Continent,y=as.numeric(Value)))+
    geom_boxplot()+
    facet_wrap(.~Stat,scales="free_y")+
    ylab("PopGen Statistic")+
    theme_bw()

ggsave("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf.pdf",
    P,
    width=16,
    height=8)

ggsave("/media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/results/Grenedalf.png",
    P,
    width=16,
    height=8)


