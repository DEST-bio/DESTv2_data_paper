# Figure 1
# panel C + PCA (fig2B) + npstat (right B)  + recomb (A) + panel D

# Libraries
library(data.table)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpattern) 
library(magick)
library(metR)
library(plyr)
library(ggthemes)
library(ggpol)
library(scales)

setwd("~/Desktop/DEST/Figure1")

##### A) Map
#### samples file
samps <- fread("DATA/dest_v2.samps_25Feb2023.csv")

samps.map <- samps %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
                               set %in% "dgn" ~ "DGN"
  ))

world <- map_data("world")

a <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    fill = "#cccccc", color = "darkgray", linewidth = 0.1
  ) +
  theme_classic(base_size = 14) +
  scale_x_longitude(breaks = seq(-160,180,40)) +
  scale_y_latitude(breaks = seq(-40,100,20),limits=c(-61,76)) +
  coord_sf() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_fill_manual(values = c("#f27f65","#8c89c1","#a5c9cc")) +
  scale_color_manual(values = c("#D05438","#615E93","#6F9FA3")) +
  labs(fill="Set", color="Set") +
  theme(legend.position = "bottom") + 
  geom_point(
    data = samps.map,
    aes(x=long,
        y=lat,
        fill = super.set, color = super.set), size = 2.5, shape = 21
  )
a

##### B) Temporal
#### Contaminantion --->
contam<-get(load("DATA/Mean.contamination.Final.Rdata"))
#### DuplicateRates --->
duprate<-get(load("DATA/DuplicateRates.all.Rdata"))
setDT(duprate)
duprate <- duprate %>%
  dplyr::select(sampleId , pcrdup)
  
pnps <-fread("DATA/pnps.predict.txt")
pnps.sub <- pnps %>%
  filter(Chrom == "genomewide") %>%
  dplyr::select(sampleId = POP, pNpS, private, Status)
  
######
pre.dat <-get(load("DATA/Miss.Cov.PCRdup.sim.joint.Rdata"))
cov.info <- pre.dat %>%
  filter(Var == "Cov") %>%
  dplyr::select(sampleId, Cov = Value)
 

Miss.info <- pre.dat %>%
  filter(Var == "Miss") %>%
  dplyr::select(sampleId, Miss = Value)
 
metadata.seq <- full_join(contam, duprate) %>%
  full_join(pnps.sub) %>%
  full_join(cov.info) %>%
  full_join(Miss.info) %>%
  full_join(samps)

metadata.seq <- metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  ))


metadata.seq.df <- metadata.seq %>%
  filter(year >= 2009)

b <- ggplot(metadata.seq.df, aes(
  y= lat,
  x= jday,
  group=city,
  color=super.set,fill=super.set
)) +
  geom_line(alpha=0.5, show.legend = F) +
  geom_point(show.legend = F,
             aes(fill = super.set, color = super.set), size = 1.5, shape = 21
) +
  scale_y_latitude(breaks = seq(-40,100,20),limits=c(-29,65),name = "Latitude") +
  scale_x_continuous(limits=c(100,350),breaks = seq(100,400,100)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9)) +
  labs(x="Julian day", y = "Latitude")+
  facet_grid(~year)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_fill_manual(values = c("#f27f65","#8c89c1","#a5c9cc")) +
  scale_color_manual(values = c("#D05438","#615E93","#6F9FA3")) +
  labs(fill="Set", color="Set") +
  theme(legend.position = "bottom")
b

##### C) PCA

dat.for.PCA <- metadata.seq %>%
  filter(set != "dgn") %>%
  mutate(collapse = case_when(collector == "Fournier-Level et al" ~ "Yes", TRUE ~ "No")) %>%
  mutate(super.set = case_when(set %in% c("DrosEU","DrosRTEC") ~ "DEST 1.0",
                               set %in% c("cville","dest_plus","DrosEU_3","DrosEU_3_sa") ~ "DEST 2.0",
  )) %>%
  group_by(super.set)

QC.pca <- dat.for.PCA %>% 
  dplyr::select(MAPPED_eff, SimCont.Norm, pcrdup, pNpS, private, Cov, Miss) %>% .[,-1] %>%
  PCA(graph = F)
  

Metadata_final <- QC.pca$ind$coord %>%
  as.data.frame() %>%
  cbind(dat.for.PCA) %>% 
  mutate(Recomendation = case_when(Status == "Exclude" ~ "Abnormal Pn/Ps",
                                   SimCont.Norm > 0.15 ~ "High contamination",
                                   Miss > 0.30  ~ "High missing data",
                                   collapse == "Yes" &  Miss < 0.30 ~ "Collapse",
                                   TRUE ~ "Pass"
  ))

Metadata_final$Recomendation <- factor(Metadata_final$Recomendation, levels=c("Pass", "Collapse", "High contamination", "Abnormal Pn/Ps", "High missing data"), labels=c("Pass", "Collapse", "High\ncontamination", "Abnormal\nPn/Ps", "High\nmissing data"))
c <-  
  ggplot(Metadata_final,
         aes(x=Dim.1,
             y=Dim.2,
             fill = Recomendation,
             shape = Recomendation,
             color = Recomendation,
  )) +
  scale_shape_manual(values = 21:25) +
  #scale_fill_gradient2(midpoint = 130, mid = "yellow", low = "blue", high = "red") +
  geom_point(size = 2.5) +
  theme_bw(base_size = 14) + labs(x="PC 1 (39.5%)", y = "PC 2 (18.6)",fill="",color="",shape="")  +
  theme(axis.text = element_text(color = "black")) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_manual(values=c("#284d42","#643827", "#384051", "#5C374E", "#425621")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  theme(axis.text=element_text(size=14), legend.text = element_text(size=8), legend.title = element_text(size=8))

c

# D) npstat

#load("DATA/pop_genome_stats_pi_m.RData")
load("DATA/diversity_dest2.Rdata")
#pop_genome_stats_pi_m$Continent<-factor(pop_genome_stats_pi_m$Continent,labels = c("Europe", "Asia", "Africa", "North\nAmerica","Oceania","South\nAmerica"))
d<-ggplot(data=divm, aes(x=continentf, y=Pi, color=continentf,fill=continentf))+
#   #geom_boxplot(data=pop_genome_stats_pi_m,aes(x=Continent,y=value,fill=Continent), show.legend = F)+
   geom_boxjitter( jitter.shape = 21,  
                   outlier.color = NA, errorbar.draw = TRUE, show.legend = F)+
  ylab(expression("Nucleotide diversity ("~paste(pi)~")"))+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        text=element_text(size=14,color="black"),axis.text = element_text(color="black")) + labs(x="")   +
  #scale_fill_gdocs() + scale_color_gdocs() 
  scale_fill_brewer(palette = "Spectral") +
  scale_color_manual(values=c("#7f2831","#a3492d","#a06e3e","#9b8a5a","#8e8e6b","#00687f","#6c8c67","#407262","#1a4863"))
d

# E) Recombination

infile <- paste0("DATA/2L.100k.txt")
data <- read.table(infile, h=T)
attach(data)

# Reshape the data from wide to long format
datalong <- gather(data, key = "variable", value = "value", -pos)

# Create the line plot
line_color <- "grey"
line_color2 <- "black"
ci_color <- "lightblue"  # Set the color for the confidence interval


chromosomes <- c("2L", "2R", "3L", "3R")
#chromosomes <- c("2L", "2R", "3L", "3R", "2L.mod", "2R.mod", "3L.mod", "3R.mod")
chromosome_lines <- list(
  "2L" = c(2225744, 13154180),
  "2R" = c(15391154, 20276334),
  "3L" = c(3173046, 16308841),
  "3R" = c(21406917, 29031297, 16432209, 24744010, 20096867, 32079331, 11750567, 26140370),
  "X" = c(13625736, 19579328, 17828912, 19593711),
  "2L.mod" = c(2225744, 13154180),  
  "2R.mod" = c(15391154, 20276334),
  "3L.mod" = c(3173046, 16308841),
  "3R.mod" = c(21406917, 29031297, 16432209, 24744010, 20096867, 32079331, 11750567, 26140370),
  "X.mod" = c(13625736, 19579328, 17828912, 19593711)
)

chromosome <- "2L"

vertical_lines <- chromosome_lines[[chromosome]]

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
e <- ggplot(datalong, aes(x = pos, y = value)) +
  annotate("rect", xmin = 2225744, xmax = 13154180, ymin = -Inf, ymax = Inf, fill="#ffe4db",alpha=0.5) +
  annotate("text", x = 6955744, y = 0.000000035, label = "In(2L)t",fontface = 'italic') +
  geom_line(color = line_color, aes(group = variable), linewidth = 0.1, na.rm=T) +
  stat_summary(fun.data = "mean_cl_normal", geom = "ribbon", alpha = 0.6, fill = ci_color, color = "transparent", na.rm=T) +
  stat_summary(fun = "mean", geom = "line", color = "#41161B", linewidth = 0.5, na.rm = TRUE) +
  labs(x = paste0("Chromosome ", chromosome), y = "rec") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black', linewidth = 0.3)) +
  theme(legend.position = "none", panel.background = element_blank(), panel.border = element_blank()) +
 # geom_vline(xintercept = vertical_lines, linetype = "dotted", color = "red", linewidth = 0.3, na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(x="Chromosome 2L", y = "Recombination rate (c/bp)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_y_continuous(label=scientific_10)
e

composite <- cowplot::plot_grid(a,b,plot_grid(c,d,ncol=2,align="tb",labels=c("C", "D")),e,labels=c("A","B","","E"),ncol=1,align="lrtb")

ggsave("figure1_v2.pdf",plot=composite, height = 15, width = 8)
