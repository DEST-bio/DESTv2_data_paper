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
  filter(pop1 != "Karensminde") %>%
  filter(pop1 != "Vyshneve")

fst.winter.between %>%
  group_by(pop1) %>%
  slice_max(ydelta,  with_ties = FALSE) %>% 
  as.data.frame()

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
pops <- unique(fst.winter.between$pop1)

over_ests_fst = 
  foreach(i = pops, .combine = "rbind")%do%{
    
    tmp <- fst.winter.between %>%
      filter(pop1 == i) 
    
    year_uniq <- unique(tmp$y1)

    foreach(k=year_uniq[-(length(year_uniq))], 
            .errorhandling = "remove",
            .combine = "rbind")%do%{
              y1.k = k
              y2.k = k+1
              message(paste(y1.k, y2.k, sep = " "))
              
              tmp %>%
                filter(y1 == y1.k) %>%
                filter(ydelta %in% 0:1) %>%
                group_by(ydelta) %>%
                summarise(mFst = mean(FST)) ->
                tmp2
              
              inner = tmp2$mFst[which(tmp2$ydelta == 0)]
              overw = tmp2$mFst[which(tmp2$ydelta == 1)]
              
              data.frame(pop=i,
                         year=k,
                         inner=inner,
                         overw=overw,
                         delta=abs(overw/inner))
            }
  }

###


###
over_ests_fst %>%
  group_by(pop) %>%
  summarise(mI = median(inner),
            mO = median(overw)) %>%
  ggplot(aes(
    x=pop,
    y=abs(mI-mO),
    fill = pop
  )) +
  geom_bar(width = 0.8, stat = "identity") +
  scale_fill_brewer(palette = "Accent") + 
  theme_bw() + coord_flip() ->
  delta.bar
ggsave(delta.bar, 
       file = "delta.bar.pdf", 
       w = 4, h = 3)


over_ests_fst %>%
  select(pop, inner, overw) %>%
  melt(id = c("pop")) %>%
  ggplot(aes(
    x=pop,
    y=value,
    color=variable
  )) +
  geom_boxplot(width = 0.5) + theme_bw() +
  scale_color_brewer(palette = "Set1")+ coord_flip() ->
  delta.box2
ggsave(delta.box2, 
       file = "delta.box2.pdf", 
       w = 4, h = 3)

####
fst.winter.between %>%
  dplyr::select(#residual.t.loc,
    #ydelta.test,
    day_diff,
    pop1,
    #T.mean,
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
  scale_color_brewer(palette = "Accent") + 
  geom_smooth(method = "lm", se =T)->
  Twinbet.ti.plot
ggsave(Twinbet.ti.plot, 
       file = "Twinbet.ti.plot.pdf", 
       w = 6, h = 3)

########
########
########

####
fst.winter.between %>%
  group_by(ydelta==0, pop1, ydelta) %>%
  filter(ydelta <= 7) %>%
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
       w = 6, h = 3)


####

test_fst_p = 
  foreach(i = pops, .combine = "rbind")%do%{
    
    tmp <- fst.winter.between %>%
      filter(pop1 == i)
    
    within <- tmp$FST[which(tmp$ydelta==0)]
    between <- tmp$FST[which(tmp$ydelta!=0)]
    wilcox.test(within, between) -> ttets
    
    data.frame(pop=i,
               pval=ttets$p.value)
  }
