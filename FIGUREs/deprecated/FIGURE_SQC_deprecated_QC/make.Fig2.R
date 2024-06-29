### MAke Figure 2.yes
### 

library(data.table)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpattern) ##--> remotes::install_github("coolbutuseless/ggpattern")
library(magick)

setwd("/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE2/")
#### samps ---> 
samps <- fread("../FIGURE1/code/dest_v2.samps_25Feb2023.csv")

####

#### Contaminantion --->
contam<-get(load("../FIGURE1/code/Mean.contamination.Final.Rdata"))
#### DuplicateRates --->
duprate<-get(load("../FIGURE1/code/DuplicateRates.all.Rdata"))
setDT(duprate)
duprate %>%
  dplyr::select(sampleId , pcrdup) ->
  duprate
#### PnPs ---->
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/misc/PnPs/results/classify_pops.txt?token=GHSAT0AAAAAAB7U52CNPDJUDB7GWR4N27GQZCAF6WA -O pnps.predict.txt")
pnps <-fread("../FIGURE1/code/pnps.predict.txt")
pnps %>%
  filter(Chrom == "genomewide") %>%
  dplyr::select(sampleId = POP, pNpS, private, Status) ->
  pnps.sub
######
pre.dat <-get(load("../FIGURE1/code/Miss.Cov.PCRdup.sim.joint.Rdata"))
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

####
####

metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) %>%
  mutate(set = case_when(set == "DrosEU_3_sa" ~ "DrosEU_3",
                         TRUE ~ set))->metadata.seq

#### CREATE FIGURE 2
####
metadata.seq  %>%
  group_by(set) %>%
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
  theme_bw() +
  ylab("") +
  xlab("Sample") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())->
  PanelA

#ggsave(QC.figure, file = "QC.figure.pdf", w = 5, h = 6)

#####

metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) %>%
  group_by(super.set) -> dat.for.PCA

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
  PanelB


### Panel C

dat.for.PCA %>% 
  dplyr::select(MAPPED_eff, SimCont.Norm, pcrdup, pNpS, private, Cov, Miss) %>% .[,-1] %>%
  PCA(graph = F) ->
  QC.pca

fviz_pca_var(QC.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) -> PanelC

###

ggsave(PanelA, file = "PanelA.pdf",
       w= 5, h = 5.5)

ggsave(PanelB, file = "PanelB.pdf",
       w= 4, h = 3.0)

ggsave(PanelC, file = "PanelC.pdf",
       w= 4, h = 3.0)
