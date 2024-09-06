#### FET enrichment across tests

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(ggvenn)
#library(GenomicRanges)
#library(rtracklayer)


#### PROCEED TO CHECKPOINT .....

#### Bergland
berg_files =
system("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/Overlap_analysis/bergland2016_LifOve",
       intern = T)

bergland2016_dm6 =
  foreach(i=berg_files,
          .combine = "rbind")%do%{
            tmp <- get(load(
              paste("/gpfs2/scratch/jcnunez/DEST2.0_analysis/Overlap_analysis/bergland2016_LifOve",
                    i, sep = "/"
              )
            ))
          }

bergland2016_dm6 %>%
  group_by(chr, pos) %>%
  arrange(chr, pos) ->
  bergland2016_dm6_srt

save(bergland2016_dm6_srt,
     file = "/netfiles/nunezlab/D_melanogaster_resources/Datasets/2016.Bergland/bergland2016_dm6_srt.Rdata")

#####
machad_files =
  system("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/Overlap_analysis/Machado2021_LifOve",
         intern = T)

Machado2021_dm6 =
  foreach(i=machad_files,
          .combine = "rbind")%do%{
            tmp <- get(load(
              paste("/gpfs2/scratch/jcnunez/DEST2.0_analysis/Overlap_analysis/Machado2021_LifOve",
                    i, sep = "/"
              )
            ))
          }

Machado2021_dm6 %>%
  group_by(chr, pos) %>%
  arrange(chr, pos) ->
  Machado2021_dm6_srt

save(Machado2021_dm6_srt,
     file = "/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.Machado/Machado2021_dm6_srt.Rdata")

##### LOAD -- CHECKPOINT

Machado2021_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.Machado/Machado2021_dm6_srt.Rdata"))
bergland2016_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2016.Bergland/bergland2016_dm6_srt.Rdata"))

c2_dat="/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest2_seasonality/dest2_glm_baypass_annotation_pod.podOutpuToo.Rdata"
load(c2_dat)

m %<>%
  mutate(SNP_id_dm6 = paste(chr,pos, sep = "_"))

Reduce(intersect, 
       list(Machado2021_dm6_srt$SNP_id_dm6,
            bergland2016_dm6_srt$SNP_id_dm6,
            m$SNP_id_dm6)
       ) -> common_SNPs

Machado2021_dm6_srt[which(Machado2021_dm6_srt$SNP_id_dm6 %in% common_SNPs),] -> machado_commons
bergland2016_dm6_srt[which(bergland2016_dm6_srt$SNP_id_dm6 %in% common_SNPs),] -> bergland_commons
m[which(m$SNP_id_dm6 %in% common_SNPs),] -> c2_commons

#### Begin comparison for FET ....
left_join(c2_commons[,c("SNP_id_dm6","C2_std_median")],
          bergland_commons[,c("SNP_id_dm6","sfsfsfX.p")],
          ) %>%
  left_join( machado_commons[,c("seas.p","SNP_id_dm6")]) %>%
  .[complete.cases(.),]->
  joint_object

##### c2 vs Berg...

one_per_c2 = thrs.ag$C2_thr[6]
b.tresh=0.3 
m.tresh=0.004

###
joint_object %>%
  filter(C2_std_median > one_per_c2 & sfsfsfX.p < b.tresh) %>% dim(.) %>% .[1] -> b.both
joint_object %>%
  filter(C2_std_median > one_per_c2 & sfsfsfX.p > b.tresh) %>% dim(.) %>% .[1] -> b.c2_only
joint_object %>%
  filter(C2_std_median < one_per_c2 & sfsfsfX.p < b.tresh) %>% dim(.) %>% .[1] -> b.berg_only
joint_object %>%
  filter(C2_std_median < one_per_c2 & sfsfsfX.p > b.tresh) %>% dim(.) %>% .[1] -> b.neither

#b.both+b.c2_only+b.berg_only+b.neither == dim(joint_object)[1]

data.frame(
  test = "vs.Bergland2016",
  scope = "genome_wide",
  both=b.both,
  c2_only=b.c2_only,
  other_only=b.berg_only,
  neither=b.neither,
  sum=b.both+b.c2_only+b.berg_only+b.neither
) -> b.all

### now per chromosome
b.chrs =
  foreach(chr.i = c("2L","2R","3R","3L"),
          .combine = "rbind")%do%{
    
    flt = joint_object[grep(chr.i, joint_object$SNP_id_dm6),]
      
          flt %>%
     filter(C2_std_median > one_per_c2 & sfsfsfX.p < b.tresh) %>% dim(.) %>% .[1] -> b.both
            flt %>%
     filter(C2_std_median > one_per_c2 & sfsfsfX.p > b.tresh) %>% dim(.) %>% .[1] -> b.c2_only
            flt %>%
     filter(C2_std_median < one_per_c2 & sfsfsfX.p < b.tresh) %>% dim(.) %>% .[1] -> b.berg_only
            flt %>%
     filter(C2_std_median < one_per_c2 & sfsfsfX.p > b.tresh) %>% dim(.) %>% .[1] -> b.neither
    
          data.frame(
              test = "vs.Bergland2016",
              scope = chr.i,
              both=b.both,
              c2_only=b.c2_only,
              other_only=b.berg_only,
              neither=b.neither,
              sum=b.both+b.c2_only+b.berg_only+b.neither
            ) 
            
          }

###### Machado

joint_object %>%
  filter(C2_std_median > one_per_c2 & seas.p < m.tresh) %>% dim(.) %>% .[1] -> m.both
joint_object %>%
  filter(C2_std_median > one_per_c2 & seas.p > m.tresh) %>% dim(.) %>% .[1] -> m.c2_only
joint_object %>%
  filter(C2_std_median < one_per_c2 & seas.p < m.tresh) %>% dim(.) %>% .[1] -> m.other_only
joint_object %>%
  filter(C2_std_median < one_per_c2 & seas.p > m.tresh) %>% dim(.) %>% .[1] -> m.neither


data.frame(
  test = "vs.Machado2021",
  scope = "genome_wide",
  both=m.both,
  c2_only=m.c2_only,
  other_only=m.other_only,
  neither=m.neither,
  sum=m.both+m.c2_only+m.other_only+m.neither
) -> m.all

m.chrs =
  foreach(chr.i = c("2L","2R","3R","3L"),
          .combine = "rbind")%do%{
            
            flt = joint_object[grep(chr.i, joint_object$SNP_id_dm6),]
            
            flt %>%
              filter(C2_std_median > one_per_c2 & seas.p < m.tresh) %>% dim(.) %>% .[1] -> m.both
            flt %>%
              filter(C2_std_median > one_per_c2 & seas.p > m.tresh) %>% dim(.) %>% .[1] -> m.c2_only
            flt %>%
              filter(C2_std_median < one_per_c2 & seas.p < m.tresh) %>% dim(.) %>% .[1] -> m.other_only
            flt %>%
              filter(C2_std_median < one_per_c2 & seas.p > m.tresh) %>% dim(.) %>% .[1] -> m.neither
            
            data.frame(
              test = "vs.Machado2021",
              scope = chr.i,
              both=m.both,
              c2_only=m.c2_only,
              other_only=m.other_only,
              neither=m.neither,
              sum=m.both+m.c2_only+m.other_only+m.neither
            ) 
            
          }

######
FET_pre_tables = rbind(b.all, b.chrs, m.all, m.chrs)

FET_output=
foreach(i=1:dim(FET_pre_tables)[1],
        .combine = "rbind")%do%{
          
          tmp <- FET_pre_tables[i,]
          
          fet_table <-
            matrix(c(tmp$both, 
                     tmp$other_only, 
                     tmp$c2_only, 
                     tmp$neither),
                   nrow = 2)
          fisher.test(fet_table) -> fet.o
          
          data.frame(
            test=tmp$test,
            scope=tmp$scope,
            both=tmp$both,
            c2_only=tmp$c2_only,
            other_only=tmp$other_only,
            neither=tmp$neither,
            fet.p = fet.o$p.value,
            fet.OR = fet.o$estimate,
            lci=fet.o$conf.int[1],
            uci=fet.o$conf.int[2]
          )
        }

save(FET_output, file = "FET_output.Rdata")

FET_output %>%
  ggplot(
    aes(
      x=scope,
      y=log2(fet.OR),
      ymin=log2(lci),
      ymax=log2(uci)
    )
  ) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 3, shape = 21, fill = "grey") +
  facet_grid(variable~test) + theme_bw() ->
  FET_output.plot
ggsave(FET_output.plot, file = "FET_output.plot.pdf",
       w = 7, h = 3)



