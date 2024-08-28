#### SNP-wise analyses -- PE Data

#### ==> Libs
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(ggvenn)

#### ===> Load PE's data
PE_snps <- fread("erickson_top1.txt") %>% mutate(SNP_id3 = paste(chr,dm3.pos, sep = "_"))

#### ===> Lift-over data
library(GenomicRanges)
library(rtracklayer)
chain_3to6="/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain"
chainObject <- import.chain(chain_3to6)

##iterate over positions
lift_over_top_snps=
foreach(i=1:dim(PE_snps)[1],
        .combine = "rbind")%do%{
          
          message(paste(i, dim(PE_snps)[1], sep = "|"))
          tmp = PE_snps[i]
          chr.i=paste("chr", tmp$chr, sep = "")
          pos.i=tmp$dm3.pos
          
          grObject <- GRanges(seqnames=c(chr.i), 
                              ranges=IRanges(start=pos.i, end=pos.i))
          results <- as.data.frame(liftOver(grObject, chainObject))
          
          data.frame(
            chr=tmp$chr,
            dm3.pos=tmp$dm3.pos,
            dm6.pos=results$start
          )
        }

save(lift_over_top_snps, file = "PE_outliers.lift_over_top_snps.Rdata")
# run liftOver
#### Load data 
load("PE_outliers.lift_over_top_snps.Rdata")
lift_over_top_snps %>%
  mutate(SNP_id3 = paste(chr,dm3.pos, sep = "_")) %>%
  mutate(SNP_id6 = paste(chr,dm6.pos, sep = "_")) %>%
  filter(chr %in% c("2L","2R","3L","3R"))  ->
  lift_over_top_snps_auts

######## C2 analyses
######## C2 analyses
######## C2 analyses
######## C2 analyses
######## C2 analyses

c2_dat="/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest2_seasonality/dest2_glm_baypass_annotation_pod.podOutpuToo.Rdata"

load(c2_dat)
m[cont.p < 0.01] -> top1_C
names(top1_C)[3] = "dm6.pos"
top1_C  %<>%
  mutate(SNP_id6 = paste(chr,dm6.pos, sep = "_"))

####
x <- list(
  C2 = unique(top1_C$SNP_id6), 
  PE = unique(lift_over_top_snps_auts$SNP_id6)
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) -> venn.plot

ggsave(venn.plot, file = "venn.plot.pdf")


##730625 2L
##583068 2R
##690750 3L
##777100 3R
##2,781,543
