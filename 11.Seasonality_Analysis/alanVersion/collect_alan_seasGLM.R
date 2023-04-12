#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions

#####
#####
#####
#####
#####
### open metadata
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

samps <- fread("./dest_v2.samps_25Feb2023.csv")
### open filtering object
#keep.remove <- get(load("/project/berglandlab/DEST2.0_working_data/keep.fail.samps.Rdata"))
#pass.filt <- filter(keep.remove, keep == "PASS")$sampleId
#### open seasonal pairs
#seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
#set.samps <- filter(seasonal.sets, sampleId %in%  pass.filt)$sampleId

#### Central files
files.seas = system("ls ./GLM_newods_ALAN_APR122023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_newods_ALAN_APR122023/"

#### binning analysis
seas.p.bin = 
  foreach(fil = files.seas, .combine = "rbind", .errorhandling = "remove" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/")))
    
    return(tmp)
    
  }

#### #### #### 
#### 
#### 
#### 
#### Collect into permutation sets
#### Collect into permutation sets
#### Collect into permutation sets
#### Collect into permutation sets

outfolder = "/project/berglandlab/jcbnunez/Shared_w_Alan/GLM_ALAN_APR122023_by_PERM_quasibinomial"
### save as independent files
  foreach(p = 1:100)%do%{
    
    message(paste(p, 
                  sep = "/"))
    
    dat.flt %>% filter(perm == p 
    ) -> dat.p
    message(p)
    save(dat.p,
         file = paste(outfolder,
                      paste("GLM_out.perm_id", p, "SeasAlanQusiBinom", "Rdata", sep = "."),
                      sep = "/"
                      ))
  }

#### Collect into permutation sets -- end
#### Collect into permutation sets -- end
#### Collect into permutation sets -- end
#### 
#### 
#### 
#### #### #### 
  


#### Analyze the p-value
seas.p.bin -> dat.flt

ot_across_perms=
foreach(p = 0:100, .combine = "rbind" )%do%{
  #foreach(chr = c("2L","2R","3L","3R"), .combine = "rbind")%do%{
    message(paste(p, 
                  #chr, 
                  sep = "/"))
    
    dat.flt %>% filter(perm == p 
                       #& chr == chr
                       ) -> dat.0
    hist(dat.0$p_lrt, breaks = 10) -> hist.dat
    data.frame(hist.dat$mids ,    hist.dat$counts) %>%
      mutate(perm = p, 
             #chr = chr
             )-> ot
    return(ot)
}
#}
#^^^ uncomment for CHR


### Plot

setDT(ot_across_perms)

ot_across_perms %>%
  ggplot(
    aes(
      x=hist.dat.mids,
      y=hist.dat.counts,
      color = perm == 0,
      group = perm
    )
  ) + 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +  geom_line(aes(alpha = perm == 0) ) +
  geom_vline(xintercept = 0.05) +
  scale_alpha_manual(values = c(0.05, 1)) ->
  pvals.lines
ggsave(pvals.lines, file = "pvals.lines.pdf", w = 9, h = 4)

  ggplot() +
  geom_boxplot(
    data = ot_across_perms[perm != 0],
    aes(x= as.factor(hist.dat.mids),
        y=hist.dat.counts),
    fill = "grey", size = 0.5, outlier.shape = NA
  ) +
    geom_point(
      data = ot_across_perms[perm == 0],
      aes(x= as.factor(hist.dat.mids),
          y=hist.dat.counts),
      size = 3,
      color = "red"
    ) + 
    #facet_wrap(~chr) +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) -> p.vals.filts
  
  ggsave(p.vals.filts, file = "p.vals.filts.pdf", w = 9, h = 4)

