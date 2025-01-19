### Plot FST analysis spac. temp
### module load Rgeospatial

library(tidyverse)
library(magrittr)
library(foreach)
library(data.table)
library(foreach)
library(rnaturalearth)

library(nasapower)
library(sp)
library(lubridate)
library(stringr)
library(car)

library(FactoMineR)
library(factoextra)
library(segmented)
library(gmodels)
library(tibble)
library(tidybulk)



#### Many winters
#### Many winters
file <- "/gpfs2/scratch/jcnunez/DEST2.0_analysis/spa_temp_fst/fst.winter.newTimePops.Rdata"
load(file)
setDT(fst.winter.between) 

fst.winter.between %<>%
  filter(pop1 != "Karensminde")


####
#anova(lm(logit(abs(FST)) ~ pop1*day_diff+T.min+T.max+T.mean, data = fst.winter.between))
#lm(logit(abs(FST)) ~ pop1*day_diff, data = fst.winter.between) -> mod.tim.loc

cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.max)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.min)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$T.mean)
cor.test(logit(abs(fst.winter.between$FST)), fst.winter.between$day_diff)

foreach(i=unique(fst.winter.between$pop1),
        .combine = "rbind")%do%{
          
          message(i)
          tmp = fst.winter.between %>%
            filter(pop1 == i)
          
          cor.test(logit(abs(tmp$FST)), tmp$T.mean)-> tp
          cor.test(logit(abs(tmp$FST)), tmp$T.min)-> tmin
          cor.test(logit(abs(tmp$FST)), tmp$T.max)-> tmax
          cor.test(logit(abs(tmp$FST)), tmp$day_diff)-> timep
          
          data.frame(
            pop=i,
            meanT.p = tp$p.value,
            minT.corr = tmin$estimate,
            minT.p = tmin$p.value,
            maxT.corr = tmax$estimate,
            maxT.p = tmax$p.value,
            meanT.corr = tp$estimate,
            Time.p = timep$p.value,
            Time.corr = timep$estimate
          )
          
        }

####
fst.winter.between %>%
  dplyr::select(#residual.t.loc,
    #ydelta.test,
    day_diff,
    pop1,
    T.mean,
    FST,
    #T.min,
    #T.max
  ) %>%
  melt(id=c("FST","pop1"
  ) ) %>% 
  ggplot(
    aes(
      x=(value),
      y=logit(abs(FST)),
      color = pop1,
    )
  ) +
  ylab(expression(paste(#"residual", 
    F[ST]))) +
  xlab("value between 2 samples") +
  geom_point(alpha = 0.1) +
  facet_grid(.~variable,
             #ncol = 2,
             scales = "free_x"
  ) +
  theme_bw() +
  geom_smooth(method = "lm", se =T)->
  Twinbet.ti.plot
ggsave(Twinbet.ti.plot, 
       file = "Twinbet.ti.plot.pdf", 
       w = 7, h = 3)

########
########
########

fst.winter.between %>%
  dplyr::select(#residual.t.loc,
    #ydelta.test,
    day_diff,
    pop1,
    T.mean,
    FST,
    #T.min,
    #T.max
  ) %>%
  melt(id=c("FST","pop1"
  ) ) %>% 
  ggplot(
    aes(
      x=(value),
      y=logit(abs(FST)),
      color = pop1,
    )
  ) +
  ylab(expression(paste(#"residual", 
    F[ST]))) +
  xlab("value between 2 samples") +
  #geom_point(alpha = 0.1) +
  facet_wrap(.~variable,
             ncol = 2,
             scales = "free_x"
  ) +
  theme_bw() +
  geom_smooth(method = "lm", se =T)->
  Twinbet.ti.plot.pop
ggsave(Twinbet.ti.plot.pop, 
       file = "Twinbet.ti.plot.pop.pdf", 
       w = 4, h = 3)

####
fst.winter.between %>%
  group_by(ydelta==0, pop1, ydelta) %>%
  filter(ydelta <= 5) %>%
  summarise(mFST = mean(FST)) %>%
  dcast(pop1+ydelta~`ydelta == 0`, value.var = "mFST") ->
  delta_df

pops <- unique(delta_df$pop1)

deltas_fst = 
foreach(i = pops, .combine = "rbind")%do%{
  
  tmp <- delta_df %>%
    filter(pop1 == i)
  
  within <- tmp$`TRUE`[complete.cases(tmp$`TRUE`)]
  between <- tmp$`FALSE`[complete.cases(tmp$`FALSE`)]
  deltas <- tmp$ydelta[-1]
    
  data.frame(pop=i,
             deltas=deltas,
             delta=between-within)
}

deltas_fst %>%
  ggplot(aes(
    y=delta,
    x=deltas,
    color=pop,
    linetype = pop=="Yesiloz"
  )) + geom_line() + theme_bw()->
  delta.line
ggsave(delta.line, 
       file = "delta.line.pdf", 
       w = 4, h = 3)
