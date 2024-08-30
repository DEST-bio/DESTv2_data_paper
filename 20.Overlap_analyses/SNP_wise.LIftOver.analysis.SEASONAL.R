#### SNP-wise analyses -- Seasonal Data

#### ==> Libs
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(ggvenn)
library(GenomicRanges)
library(rtracklayer)
chain_3to6="/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain"
chainObject <- import.chain(chain_3to6)


######## C2 analyses

c2_dat="/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest2_seasonality/dest2_glm_baypass_annotation_pod.podOutpuToo.Rdata"

load(c2_dat)

trehsh=thrs.ag$C2_thr[6]

m %>%
  filter(C2_std_median >trehsh ) %>%
  mutate(SNP_id_dm6 = paste(chr,pos, sep = "_") ) ->
  m.c2.outlrs

#### Bergland data
load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2016.Bergland/6d_data.Rdata")
### note in DM3
  p %>%
  filter(sfsfsfX.q<.3) ->
    p_seasonal_hits

  lift_over_p2016_seas=
    foreach(i=1:dim(p_seasonal_hits)[1],
            .combine = "rbind",
            .errorhandling = "remove")%do%{
              
              message(paste(i, dim(p_seasonal_hits)[1], sep = "|"))
              tmp = p_seasonal_hits[i,]
              
              chr.i=paste("chr", tmp$chr, sep = "")
              pos.i=tmp$pos
              
              grObject <- GRanges(seqnames=c(chr.i), 
                                  ranges=IRanges(start=pos.i, end=pos.i))
              results <- as.data.frame(liftOver(grObject, chainObject))
              
              data.frame(
                chr=tmp$chr,
                pos=tmp$pos,
                dm6.pos=results$start
              )
            }
  
  lift_over_p2016_seas %>%
    mutate(SNP_id_dm6 = paste(chr,dm6.pos, sep = "_") ) %>%
    right_join(p_seasonal_hits) -> 
    lift_over_p2016_seas.LOft
  
  save(lift_over_p2016_seas.LOft, 
       file = "lift_over_p2016_seas.LOft.Rdata")

  ######  
  # Machado
  mach_dat <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.Machado/seas_glm_NOTswitch_clean.txt")
  mach_dat %>%
    filter(seas.p < 0.004) ->
    machado_seas_snps
  ### 0.004 is the top 1% as per the paper
  
  lift_over_machado_seas_snps=
    foreach(i=1:dim(machado_seas_snps)[1],
            .combine = "rbind",
            .errorhandling = "remove")%do%{
              
              message(paste(i, dim(machado_seas_snps)[1], sep = "|"))
              tmp = machado_seas_snps[i,]
              
              chr.i=paste("chr", tmp$chr, sep = "")
              pos.i=tmp$pos
              
              grObject <- GRanges(seqnames=c(chr.i), 
                                  ranges=IRanges(start=pos.i, end=pos.i))
              results <- as.data.frame(liftOver(grObject, chainObject))
              
              data.frame(
                chr=tmp$chr,
                pos=tmp$pos,
                dm6.pos=results$start
              )
            }
  
  save(lift_over_machado_seas_snps, 
       file = "lift_over_machado_seas_snps.Rdata")
  

  ######  
  # Nunez
  
nunez <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Revision_Best_Models/temp.max;2;5.Cville.v2.glmRNP.Rdata"))  
  
  nunez %>%
    filter(perm == 0) %>%
    filter(rnp < 0.01 ) %>%
    mutate(SNP_id_dm6 = paste(chr,pos, sep = "_") ) ->
    nunez.outlrs
  
  
  ###### LOAD DATS
  bergland<-get(load("lift_over_p2016_seas.LOft.Rdata"))
  machado<-get(load("lift_over_machado_seas_snps.Rdata")) %>%
    mutate(SNP_id_dm6 = paste(chr, dm6.pos, sep = "_"))
  nunez <- nunez.outlrs
  c2_out <- m.c2.outlrs
    
  
  ####
  x <- list(
    bergland = unique(bergland$SNP_id_dm6), 
    machado = unique(machado$SNP_id_dm6),
    nunez.o = unique(nunez.outlrs$SNP_id_dm6),
    c2_out = unique(c2_out$SNP_id_dm6)
  )

  ggvenn(
    x, 
    #fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 0.5, set_name_size = 4
  ) -> venn.plot
  
  ggsave(venn.plot, file = "venn.plot.pdf")
  

