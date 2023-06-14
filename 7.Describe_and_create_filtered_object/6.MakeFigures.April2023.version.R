
library(data.table)
library(tidyverse)
library(FactoMineR)
library(factoextra)

#### samps ---> 
samps <- fread("dest_v2.samps_25Feb2023.csv")

####

#### Contaminantion --->
contam<-get(load("/scratch/yey2sn/DEST2_analysis/filtering/Mean.contamination.Final.Rdata"))
#### DuplicateRates --->
duprate<-get(load("/scratch/yey2sn/DEST2_analysis/filtering/DuplicateRates.all.Rdata"))
setDT(duprate)
duprate %>%
  dplyr::select(sampleId , pcrdup) ->
  duprate
#### PnPs ---->
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/misc/PnPs/results/classify_pops.txt?token=GHSAT0AAAAAAB7U52CNPDJUDB7GWR4N27GQZCAF6WA -O pnps.predict.txt")
pnps <-fread("pnps.predict.txt")
pnps %>%
  filter(Chrom == "genomewide") %>%
  dplyr::select(sampleId = POP, pNpS, private, Status) ->
  pnps.sub
######
pre.dat <-get(load("Miss.Cov.PCRdup.sim.joint.Rdata"))
pre.dat %>%
  filter(Var == "Cov") %>%
  dplyr::select(sampleId, Cov = Value) -> 
  cov.info

pre.dat %>%
  filter(Var == "Miss") %>%
  dplyr::select(sampleId, Miss = Value) -> 
  Miss.info

full_join(contam, duprate) %>%
  full_join(pnps.sub) %>%
  full_join(cov.info) %>%
  full_join(Miss.info) %>%
  full_join(samps) -> 
  metadata.seq
  
metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) ->metadata.seq

save(metadata.seq, file = "fig1.fig2.Rdata")
load("fig1.fig2.Rdata")

#### CREATE FIGURE 1
world <- map_data("world")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  ) + theme_classic() -> base_world

ggsave(base_world, file = "base_world.pdf", w = 6, h = 3.5)

#### Panel of time

samps %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
                               set %in% "dgn" ~ "DGN"
  )) ->samps.map

base_world + 
  geom_point(
    data = samps.map,
    aes(x=long,
        y=lat,
        fill = super.set), size = 1.5, shape = 21
  ) +
  scale_fill_brewer(palette = "Set1") -> DEST_world

ggsave(DEST_world, file = "DEST_world.pdf", w = 7, h = 3.5)

###
metadata.seq %>%
  filter(year >= 2009) %>%
  ggplot(aes(
    y= lat,
    x= jday,
    group=city,
    color=super.set
  )) +
  geom_line() +
  geom_point(size = 1.2) +
  ylim(20,65) +
  xlim(100,350) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
  facet_grid(~year) ->
  time.plot
ggsave(time.plot, file = "time.plot.pdf", w = 8.5, h = 3.5)


## UpSet plot

samps %>%
  group_by(year, city) %>%
  filter(set != "dest_plus") %>% 
  summarize(N = n()) %>%
  .[complete.cases(.),] %>%
  filter(year > 2006) %>%
  filter(N >= 2) %>%
  ggplot(aes(
    y=city,
    x=year,
    fill = N
  )) + geom_point(shape = 22) +
  scale_fill_gradient2(midpoint = 7.5, low = "blue", high = "red", mid = "grey") +
  theme_classic() ->
  time_resolution
ggsave(time_resolution, file = "time_resolution.pdf", h = 5, w = 3.5)



#### CREATE FIGURE 2
####
metadata.seq  %>%
  group_by(super.set) %>%
  mutate(rank.i = rank(Cov)) %>%
  dplyr::select(sampleId, Contamination=SimCont.Norm, Duplication=pcrdup, "Missing Data"=Miss, Coverage=Cov, Status, set, continent, super.set, "Flies Pooled"=nFlies, collapse, rank.i) %>%
  melt(id = c("sampleId","set","Status","continent","super.set", "collapse" , "rank.i")) ->
  plot.data
setDT(plot.data)

ggplot() + 
  geom_point(
    data=plot.data,
    aes(
    x=rank.i,
    y=(value),
    #color = set,
    shape = set
  ), size = 2, alpha = 0.5) +
  #geom_point(data=plot.data[collapse == "Yes"],
  #           aes(
  #             x=rank.i,
  #             y=(value),
  #             #shape = collapse,
  #             #fill = set,
  #           ), shape = 21, size = 2, alpha = 0.3, color = "black") +
  facet_grid(variable~super.set, scales = "free") +
  theme_bw() +
  ylab("") +
  xlab("Sample") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) ->
  QC.figure

ggsave(QC.figure, file = "QC.figure.pdf", w = 5, h = 6)

#####

metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) %>%
  group_by(super.set) -> dat.for.PCA

dat.for.PCA %>% 
  dplyr::select(MAPPED_eff, SimCont.Norm, pcrdup, pNpS, private, Cov, Miss) %>% .[,-1] %>%
  PCA(graph = F) ->
  QC.pca

fviz_pca_var(QC.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) -> var.pca
ggsave(var.pca, file = "var.pca.pdf", w = 6, h = 4)


QC.pca$ind$coord %>%
  as.data.frame() %>%
  cbind(dat.for.PCA) %>% 
  mutate(Recomendation = case_when(Status == "Exclude" ~ "Abnormal PNPS",
                            SimCont.Norm > 0.15 ~ "High Contamination",
                            Miss > 0.30  ~ "High Missing Data",
                            collapse == "Yes" &  Miss < 0.30 ~ "Collapse",
                            TRUE ~ "Pass"
                            )) ->
  Metadata_final

Metadata_final[,c("sampleId","SimCont.Norm","pcrdup","pNpS","private","Cov","Miss","Recomendation")] -> Recomendations

save(Recomendations, file = "Recomendations.Rdata")
fwrite(Recomendations, file = "QC.recomendations.csv")

Metadata_final %>%
  ggplot(aes(x=Dim.1,
             y=Dim.2,
             fill = Recomendation,
             shape = Recomendation,
             fill = Cov
             )) +
  scale_shape_manual(values = 21:25) +
  #scale_fill_gradient2(midpoint = 130, mid = "yellow", low = "blue", high = "red") +
  geom_point(size = 2.5) ->
  QC.pca.plot

ggsave(QC.pca.plot, file = "QC.pca.plot.pdf", w = 6, h = 4)

