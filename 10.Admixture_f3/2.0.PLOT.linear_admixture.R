# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  library(reshape2)
  
### open GDS
aus <- system("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/f3_revist/linear_admix/Australia/*", intern = T)
nam <- system("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/f3_revist/linear_admix/N.America/*", intern = T)
sam <- system("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/f3_revist/linear_admix/S.America/*", intern = T)

aus.samps = 
  foreach(fil = aus, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

nam.samps = 
  foreach(fil = nam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))

sam.samps = 
  foreach(fil = sam, .combine = "rbind")%do%{
    tmp <- get(load(fil))
  } %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate))


rbind(
aus.samps,
nam.samps,
sam.samps) ->
  linear.admix.dat.filters

save(linear.admix.dat.filters, file = "linear.admix.dat.Rdata")
load("linear.admix.dat.Rdata")
####
regression.coeffs = 
foreach(filter.i = unique(linear.admix.dat.filters$filter),
        .combine = "rbind")%do%{
        
linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "Australia" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary() -> o1
          
linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "S.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary()  -> o2

linear.admix.dat.filters %>%
  filter(filter == filter.i) %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  lm(mean.est ~ lat , data = .) %>% summary()  -> o3

rbind(
data.frame(
mod = "Australia",
lat = o1$coefficients[2,1],
p = o1$coefficients[2,4],
filter =  filter.i),
data.frame(
mod = "S.America",
lat = o2$coefficients[2,1],
p = o2$coefficients[2,4],
filter =  filter.i),
data.frame(
mod = "N.America",
lat = o3$coefficients[2,1],
p = o3$coefficients[2,4],
filter =  filter.i))

        }

regression.coeffs %>%
  arrange(mod) %>%
  mutate(signif = case_when(p < 0.05 ~ "SIG", 
                            TRUE ~ "NS"))

regression.coeffs %>%
  reshape2::dcast(filter~mod, value.var = "p")


linear.admix.dat.filters %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  arrange(lat)
####
linear.admix.dat.filters %>%
  filter(source_pop == "AFRICA") %>%
  filter(filter %in% c("noINV","All")) %>%
ggplot(aes(
  x=lat,
  y=mean.est,
  ymin = mean.est - sd.est,
  ymax = mean.est + sd.est,
  color = source_pop
)) + 
  geom_smooth(method = "lm",aes(linetype = filter)) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2.1, 
             fill = "grey",
             color = "black", aes(shape = filter)) +
  theme_bw() +
  scale_shape_manual(values = 21:22) +
  facet_grid(.~admix.set, scales = "free_x")->
  plot.admix.flt

ggsave(plot.admix.flt, file = "plot.admix.flt.pdf", w = 10, h  =3)

#### REVISION Feb 2025
linear.admix.dat.filters %>%
  filter(admix.set == "N.America" & 
           source_pop == "AFRICA") %>%
  filter(filter == "noINV")->
  NAM.DAT

L = dim(NAM.DAT)[1]

ran_samp_test =
foreach(i=1:100, .combine = "rbind")%do%{
  
  if(i==1){
    b=summary(lm(mean.est~lat, data = NAM.DAT))$coeff[2,1]
    p=summary(lm(mean.est~lat, data = NAM.DAT))$coeff[2,4]
  }
  if(i!=1){
    ran = sample(L, 50)
    b=summary(lm(mean.est~lat, data = NAM.DAT[ran,]))$coeff[2,1]
    p=summary(lm(mean.est~lat, data = NAM.DAT[ran,]))$coeff[2,4]
  }
  data.frame(
    i=i,
    b=b,
    p=p
  )
}
ran_samp_test %<>%
  mutate(test = "randomization_50")

ggplot() +
  geom_violin(
    data = filter(ran_samp_test, i != 1),
    aes(
    x=test,
    y=b
  )) + geom_point(
    data = filter(ran_samp_test, i == 1),
    aes(
      x=test,
      y=b
    ), size = 5) ->
  plot_ran_50

ggsave(plot_ran_50, file = "plot_ran_50.pdf",
       w=3, h =3)


ggplot(data = ran_samp_test,
       aes(
         y=-log10(p),
         x=b,
         color = i==1
       )) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.01)) ->
  plot_ran_50_p

ggsave(plot_ran_50_p, file = "plot_ran_50_p.pdf",
       w=4, h =3)


####
linear.admix.dat.filters %>%
  filter(source_pop == "AFRICA") %>%
  filter(filter %in% c("noINV","silent")) ->
  admix.dat.sets

adm5 = get(load("nclust.5.sampleId.cluster.Rdata"))
names(adm5)[2] = "adm5"
adm8 = get(load("nclust.8.sampleId.cluster.Rdata"))
names(adm8)[2] = "adm8"

admix.dat.sets %>%
  filter(admix.set == "N.America") %>%
  left_join(adm5) %>%
  left_join(adm8) ->
  admix.dat.sets.adms

#### Some summaries
admix.dat.sets.adms %>%
  group_by(adm8) %>%
  summarize(meanBeta = mean(mean.est))

### Plots
admix.dat.sets.adms %>%
  filter(!is.na(adm5)) %>%
  ggplot(aes(
    x=as.factor(adm5),
    y=mean.est,
    fill = as.factor(adm5)
  )) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_jitter(alpha = 0.3, color = "grey") +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~as.factor(filter)) ->
  adm5box
ggsave(adm5box, file = "adm5box.pdf", w = 4, h  =2.0)

admix.dat.sets.adms %>%
  filter(!is.na(adm8)) %>%
  ggplot(aes(
    x=as.factor(adm8) ,
    y=mean.est,
    fill = as.factor(adm8)
  )) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_jitter(alpha = 0.3, color = "grey") +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~as.factor(filter)) ->
  adm8box
ggsave(adm8box, file = "adm8box.pdf", w = 4, h  =2.0)

#######
####### pt 2. longitude analtsis
regression.coeffs.long = 
  foreach(filter.i = unique(linear.admix.dat.filters$filter),
          .combine = "rbind")%do%{
            
            linear.admix.dat.filters %>%
              filter(filter == filter.i) %>%
              filter(admix.set == "Australia" & source_pop == "AFRICA") %>%
              lm(mean.est ~ long , data = .) %>% summary() -> o1
            
            linear.admix.dat.filters %>%
              filter(filter == filter.i) %>%
              filter(admix.set == "S.America" & source_pop == "AFRICA") %>%
              lm(mean.est ~ long , data = .) %>% summary()  -> o2
            
            linear.admix.dat.filters %>%
              filter(filter == filter.i) %>%
              filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
              lm(mean.est ~ long , data = .) %>% summary()  -> o3
            
            rbind(
              data.frame(
                mod = "Australia",
                long = o1$coefficients[2,1],
                p = o1$coefficients[2,4],
                filter =  filter.i),
              data.frame(
                mod = "S.America",
                long = o2$coefficients[2,1],
                p = o2$coefficients[2,4],
                filter =  filter.i),
              data.frame(
                mod = "N.America",
                long = o3$coefficients[2,1],
                p = o3$coefficients[2,4],
                filter =  filter.i))
            
          }

regression.coeffs.long %>%
  reshape2::dcast(filter~mod, value.var = "p")

#### Plot
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

########
linear.admix.dat.filters %>%
  filter(filter == "noINV") %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA") %>%
  group_by(lat,long) %>%
  summarize(m.est = mean(mean.est)) -> NAMAncestry

## world plot
ggplot(data = world) +
  geom_sf(fill= "grey70") +
  coord_sf(xlim = c(-50, -135), ylim = c(15,50), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = NAMAncestry , 
              aes(
                x=long,
                y=lat,
                fill = m.est,
              ), alpha = 0.9, size = 3.5, shape = 21, color = "white") + 
  scale_fill_viridis_c()  ->
  WordAFR.admix

ggsave( WordAFR.admix, file = "WordAFR.admix.pdf", w = 9, h  =3)

###
NAMAncestry %>%
  ggplot(aes(
    x=long,
    y=m.est,
  )) + geom_smooth() + xlim(-50, -135) ->
  AFRhist
ggsave( AFRhist, file = "AFRhist.pdf", w = 9, h  =3)

  
