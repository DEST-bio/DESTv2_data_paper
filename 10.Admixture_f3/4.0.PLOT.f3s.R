# ijob -c 4 --mem=20G -p standard -A berglandlab_standard

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(vroom)
  library(reshape2)
  
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggExtra)
  library(foreach)
  
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"
samps <- fread(meta_git)
setDT(samps)

samps %>%
dplyr::select(focal=sampleId, focal_contnent = continent) ->
samps.focal

### plot a world  
  world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

### Load file

f3.admix.dat <- get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/f3_data/f3_meta.joint.Rdata"))

f3.admix.dat %<>%
filter(!is.na(`Z-score`)) %>%
left_join(samps.focal)

#### Begin filtering
#### 
f3.admix.dat %<>%
mutate(admix_evidence = case_when(`Z-score` < -1.65 ~ "admix",
`Z-score` >= -1.65 ~ "notadmix"))




####
####
f3.admix.dat %>%
group_by(admix_evidence, african_parent) %>%
summarize(N = n()) %>%
dcast(african_parent~admix_evidence, 
value.var = "N", fill = 0) %>%
mutate(tot = admix + notadmix) %>%
mutate(perc_admix = admix/tot*100 ) %>%
mutate(sampleId = african_parent) %>%
left_join(samps) ->
african_donor

cor.test(~ perc_admix + Cov, data = african_donor)
cor.test(~ perc_admix + nFlies, data = african_donor)



african_donor %>%
dplyr::select(perc_admix,lat,long, focal_contnent) %>%
melt(id = c("focal_contnent", "perc_admix")) %>%
ggplot(aes(
x= value, y = perc_admix)) +
geom_point() + geom_smooth(method = "lm") +
facet_grid(focal_contnent~variable) -> latLong_AF

#regressions
ggsave( latLong_AF, file = "latLong_AF.pdf", w = 9, h  =3)

## world plot
  ggplot(data = world) +
  geom_sf(fill= "grey70") +
    coord_sf(xlim = c(-25.50, 56.00), ylim = c(38.00, -39.00), expand = FALSE) + 
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.2), 
        panel.background = element_rect(fill = "white")) +
geom_jitter(data = african_donor , 
           aes(
  x=long,
  y=lat,
  fill = perc_admix,
  size = perc_admix,
), alpha = 0.9, size = 3.5, shape = 21, color = "white") + 
#+ facet_grid(~focal_contnent)
  scale_fill_viridis_c()  ->
  AFR.plot.f3admix

ggsave( AFR.plot.f3admix, file = "AFR.plot.f3admix.pdf", w = 9, h  =3)

####
####
####
#### --> EU donor
####
####

f3.admix.dat %>%
group_by(admix_evidence, european_parent) %>%
summarize(N = n()) %>%
dcast(european_parent~admix_evidence, 
value.var = "N", fill = 0) %>%
mutate(tot = admix + notadmix) %>%
mutate(perc_admix = admix/tot*100 ) %>%
mutate(sampleId = european_parent) %>%
left_join(samps) ->
eu_donor

eu_donor %>% 
group_by(country) %>%
summarize(mean.perc = mean(perc_admix)) %>%
arrange(mean.perc) %>% as.data.frame

eu_donor %>% arrange(perc_admix) %>% tail

eu_donor$perc_admix %>% mean

cor.test(~ perc_admix + lat, data = eu_donor)
cor.test(~ perc_admix + long, data = eu_donor)
cor.test(~  Cov + perc_admix, data = eu_donor)
cor.test(~ perc_admix + nFlies, data = eu_donor)
cor.test(~ Cov + lat, data = eu_donor)
cor.test(~ Cov + long, data = eu_donor)

eu_donor %>%
ggplot(aes(
x=Cov, y = perc_admix)) +
geom_point() + geom_smooth(method = "lm")  -> COV_EU
ggsave( COV_EU, file = "COV_EU.pdf", w = 9, h  =3)

### Admixed samples
f3.admix.dat %>%
group_by(admix_evidence, focal) %>%
summarize(N = n()) %>%
dcast(focal~admix_evidence, 
value.var = "N", fill = 0) %>%  
mutate(tot = admix + notadmix) %>%
mutate(perc_admix = admix/tot*100 ) %>%
mutate(sampleId = focal) %>%
left_join(samps) ->
focal_samps_f3

focal_samps_f3 %>%
group_by(continent) %>%
summarize(mean.pf3 = mean(perc_admix))

#### Bring in the linear admixture
linear_admix <-  get(load("linear.admix.dat.Rdata"))

linear_admix %>%
filter(source_pop == "AFRICA" & filter == "noINV") ->
afr.linear.admix

## joint
left_join(focal_samps_f3, afr.linear.admix) ->
admix.f3.beta

admix.f3.beta %>%
ggplot(aes(x= log10(perc_admix+0.1), y = mean.est)) +
geom_point() ->
f3.admixbeta.plot

ggsave( f3.admixbeta.plot, file = "f3.admixbeta.plot.pdf", w = 9, h  =3)

