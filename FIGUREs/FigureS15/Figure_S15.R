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
library(ggpubr)

setwd("~/Desktop/DEST/FigureS7")

sample_list <- fread("DATA/out.baypass.detpop", header = F)
colnames(sample_list) <- "sample"

samps <- fread("DATA/dest_v2.samps_8Jun2023.csv", sep = ",")
samps_subset <- samps[samps$sampleId%in%sample_list$sample]
samps_subset <- samps_subset[,c("sampleId", "country", "continent", 
                                "lat", "long")]
samps_subset <- samps_subset[match(sample_list$sample, samps_subset$sampleId),]

pca_result_2L_1_2_3 <- pca_result_2L$pop.loadings[,c(1,2,3)]
pca_result_2L_1_2_3 <- as.data.frame(pca_result_2L_1_2_3)
colnames(pca_result_2L_1_2_3) <- c("PC1", "PC2", "PC3")
pca_result_2L_1_2_3$sampleId <- rownames(pca_result_2L_1_2_3)

pca_result_df_2L <- merge(samps_subset, pca_result_2L_1_2_3, by = "sampleId")
pca_result_df_2L$type <- "2L"

pca_result_2R_1_2_3 <- pca_result_2R$pop.loadings[,c(1,2,3)]
pca_result_2R_1_2_3 <- as.data.frame(pca_result_2R_1_2_3)
colnames(pca_result_2R_1_2_3) <- c("PC1", "PC2", "PC3")
pca_result_2R_1_2_3$sampleId <- rownames(pca_result_2R_1_2_3)

pca_result_df_2R <- merge(samps_subset, pca_result_2R_1_2_3, by = "sampleId")
pca_result_df_2R$type <- "2R"


pca_result_3L_1_2_3 <- pca_result_3L$pop.loadings[,c(1,2,3)]
pca_result_3L_1_2_3 <- as.data.frame(pca_result_3L_1_2_3)
colnames(pca_result_3L_1_2_3) <- c("PC1", "PC2", "PC3")
pca_result_3L_1_2_3$sampleId <- rownames(pca_result_3L_1_2_3)

pca_result_df_3L <- merge(samps_subset, pca_result_3L_1_2_3, by = "sampleId")
pca_result_df_3L$type <- "3L"


pca_result_3R_1_2_3 <- pca_result_3R$pop.loadings[,c(1,2,3)]
pca_result_3R_1_2_3 <- as.data.frame(pca_result_3R_1_2_3)
colnames(pca_result_3R_1_2_3) <- c("PC1", "PC2", "PC3")
pca_result_3R_1_2_3$sampleId <- rownames(pca_result_3R_1_2_3)

pca_result_df_3R <- merge(samps_subset, pca_result_3R_1_2_3, by = "sampleId")
pca_result_df_3R$type <- "3R"


pca_result_aut_1_2_3 <- pca_result_aut$pop.loadings[,c(1,2,3)]
pca_result_aut_1_2_3 <- as.data.frame(pca_result_aut_1_2_3)
colnames(pca_result_aut_1_2_3) <- c("PC1", "PC2", "PC3")
pca_result_aut_1_2_3$sampleId <- rownames(pca_result_aut_1_2_3)

pca_result_df_aut <- merge(samps_subset, pca_result_aut_1_2_3, by = "sampleId")
pca_result_df_aut$type <- "Autosomes"

pca_result_df <- rbind(pca_result_df_2L, pca_result_df_2R, pca_result_df_3L, pca_result_df_3R, pca_result_df_aut)

figS7A<-ggplot(pca_result_df, aes(PC1, PC2, color=continent)) + geom_point() + facet_grid(.~type) +
  geom_point(size = 0.8) +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Spectral", name="Continent") +
  theme(legend.position = "bottom")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) 
figS7A
figS7B<-ggplot(pca_result_df, aes(PC1, PC3, color=continent)) + geom_point() + facet_grid(.~type) +
  geom_point(size = 0.8) +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Spectral", name="Continent") +
  theme(legend.position = "top")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) 
figS7B


library(tidyverse)
library(foreach)
library(reshape2)

pca_result_df$type %>% table

cors.pcs.chrs=
  foreach( ch = c("2L", "2R", "3L", "3R"), 
           .combine = "rbind")%do%{
             
             message(ch)
             pca_result_df %>%
               filter(type == ch) ->
               tmp
             
             foreach(cont = c("Europe", "North_America"), 
                     .combine = "rbind")%do%{
                       message(cont)
                       
                       tmp %>%
                         filter(continent == cont) ->
                         tmp2
                       
                       la1= cor.test(tmp2$PC1, tmp2$lat)
                       la2= cor.test(tmp2$PC2, tmp2$lat)
                       la3= cor.test(tmp2$PC3, tmp2$lat)
                       
                       lo1= cor.test(tmp2$PC1, tmp2$long)
                       lo2= cor.test(tmp2$PC2, tmp2$long)
                       lo3= cor.test(tmp2$PC3, tmp2$long)
                       
                       data.frame(
                         ch,
                         cont,
                         la1_e = la1$estimate,
                         la1_lci = la1$conf.int[1],
                         la1_uci = la1$conf.int[2],
                         
                         la2_e = la2$estimate,
                         la2_lci = la2$conf.int[1],
                         la2_uci = la2$conf.int[2],
                         
                         la3_e = la3$estimate,
                         la3_lci = la3$conf.int[1],
                         la3_uci = la3$conf.int[2],
                         
                         lo1_e = lo1$estimate,
                         lo1_lci = lo1$conf.int[1],
                         lo1_uci = lo1$conf.int[2],
                         
                         lo2_e = lo2$estimate,
                         lo2_lci = lo2$conf.int[1],
                         lo2_uci = lo2$conf.int[2],
                         
                         lo3_e = lo3$estimate,
                         lo3_lci = lo3$conf.int[1],
                         lo3_uci = lo3$conf.int[2]
                       )
                       
                     }}

cors.pcs.chrs.df <- cors.pcs.chrs %>%
  melt(id = c("ch", "cont")) %>%
  separate(variable, into = c("var", "est"), sep = "_") %>% 
  dcast(ch+cont+var~est, value.var = "value") %>%
  separate(var, into = c("vat","pc"), 2)
cors.pcs.chrs.df$vat <- factor(cors.pcs.chrs.df$vat, labels=c("Latitude", "Longitude"))
cors.pcs.chrs.df$cont <- factor(cors.pcs.chrs.df$cont, labels=c("Europe", "North America"))

figS7C <-  ggplot(cors.pcs.chrs.df,aes(
    x=vat,
    color=pc,
    y=e,
    ymin=lci,
    ymax=uci,
  )) +
  geom_hline(yintercept =  0,) +
  geom_errorbar(width = 0.1,position=position_dodge(width=1) ) + 
  geom_point(position=position_dodge(width=1)) + 
  facet_grid(cont~ch) +
    labs(x="", y="Correlation",color="Principal component")  + theme_bw(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black")) +
    theme(legend.position = "bottom")

figS7top<-ggarrange(figS7A,figS7B,nrow = 2,labels="AUTO",common.legend = T,align="hv",legend = "bottom")
figS7<-ggarrange(figS7top,figS7C,nrow=2,heights = c(1,0.7),labels=c("","C"))
ggsave(filename = "FigureS7.pdf", figS7, height = 6, width = 8)

######################
# Canonical correlation
random_pca_aut <- as.data.frame(pca_result_df_aut[,c("PC1","PC2")])
rownames(random_pca_aut) <- pca_result_df_aut$sampleId
dist_random_pca_aut <- dist(random_pca_aut)
dist_random_pca_aut_m<-as.matrix(dist_random_pca_aut)
dist_random_pca_aut_m <- as.data.frame(dist_random_pca_aut_m)

all.dat.pca.aut<-PCA.results.df[PCA.results.df$case=="all",c("sampleId","Dim.1","Dim.2")]
rownames(all.dat.pca.aut) <- all.dat.pca.aut$sampleId
all.dat.pca.aut$sampleId <- NULL
colnames(all.dat.pca.aut) <- c("PC1", "PC2")

order_samples <- rownames(dist_random_pca_aut_m)
all.dat.pca.aut <- all.dat.pca.aut[order_samples, ]
dist_all.dat.pca.aut <- dist(all.dat.pca.aut)
dist_all.dat.pca.aut_m<-as.matrix(dist_all.dat.pca.aut)
dist_all.dat.pca.aut_m <- as.data.frame(dist_all.dat.pca.aut_m)
# Convert the symmetric matrix into a long data frame
convert_to_long <- function(df) {
  # Add row names as a column
  df$row <- row.names(df)
  
  # Convert to long format
  long_df <- tidyr::pivot_longer(df, cols = -row, names_to = "column", values_to = "value")
  
  # Sort row and column to remove redundancy
  long_df <- long_df %>%
    mutate(
      row = pmin(row, column),
      column = pmax(row, column)
    ) %>%
    distinct(row, column, .keep_all = TRUE)
  
  long_df <- long_df[long_df$row!=long_df$column,]
  
  # Return the long dataframe
  return(long_df)
}

dist_random_pca_aut_m_long <- convert_to_long(dist_random_pca_aut_m)

dist_all.dat.pca.aut_m_long <- convert_to_long(dist_all.dat.pca.aut_m)

cancor(dist_random_pca_aut_m_long$value,dist_all.dat.pca.aut_m_long$value)

