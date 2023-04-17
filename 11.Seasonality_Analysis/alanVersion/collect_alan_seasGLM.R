#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions


log10_ceiling <- function(x) {
  10^(ceiling(log10(x)))
}


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
files.seas = system("ls ./GLM_omnibus_ALAN_APR122023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_omnibus_ALAN_APR122023/"

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

outfolder = "/project/berglandlab/jcbnunez/Shared_w_Alan/GLM_omnibus_ALAN_APR122023_by_PERM"
### save as independent files
  foreach(p = 0:20)%do%{
    
    message(paste(p, 
                  sep = "/"))
    
    seas.p.bin %>% filter(perm == p 
    ) -> dat.p
    save(dat.p,
         file = paste(outfolder,
                      paste("GLM_out.perm_id", p, "Omni", "Rdata", sep = "."),
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
  foreach(mfet = "PhyloRan_LocRan", .combine = "rbind"  )%do%{
    
  foreach(p = 0:20, .combine = "rbind" )%do%{
    #foreach(chr = c("2L","2R","3L","3R"), .combine = "rbind")%do%{
    message(paste(p, 
                  #chr, 
                  sep = "/"))
    
    dat.flt %>% 
      filter(perm == p 
                       & model_features == mfet
    ) -> dat.low
    
    hist(dat.low$p_lrt, breaks = 100) -> hist.dat.low
    data.frame(hist.dat.low$mids ,    hist.dat.low$counts) %>%
      mutate(perm = p, 
             model_features = mfet
      )-> ot

    return(ot)
  } }

### Plot

setDT(ot_across_perms)

ot_across_perms %>%
  ggplot(
    aes(
      x=hist.dat.low.mids,
      y=(hist.dat.low.counts),
      color = perm == 0,
      group = perm
    )
  ) + 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +  geom_line(aes(alpha = perm == 0) ) +
  geom_vline(xintercept = 0.05) +
  facet_wrap(~model_features, scales = "free_y") + 
  theme_bw() +
  xlab("Seasonal P-value") +
  ylab("Count") +
  scale_alpha_manual(values = c(0.1, 1)) ->
  pvals.lines
#ggsave(pvals.lines, file = "pvals.lines.pdf", w = 9, h = 9)
ggsave(pvals.lines, file = "pvals.lines.png", w = 5, h = 4)

####
dat.flt %>% 
  group_by(perm, model_features) %>%
  mutate(bin.id = case_when(
    p_lrt < 1e-5 ~ "e-5",
    p_lrt < 1e-4 ~ "e-4",
    p_lrt < 1e-3 ~ "e-3",
    p_lrt < 1e-2 ~ "e-2",
    p_lrt < 1e-1 ~ "e-1"
  )) ->
  dat.flt.bin.id

dat.flt.bin.id %>% 
  group_by(perm, bin.id, chr, model_features) %>%
  summarize(N = n()) %>%
  .[complete.cases(.),] -> cum.dat
setDT(cum.dat)

  ggplot() +
    geom_boxplot(
      data = cum.dat[perm != 0][model_features == c("PhyloRan_LocRan")],
    aes(x=chr, color =chr ,y=(N)),
    outlier.shape = NA
    ) +
  geom_point(
    data = cum.dat[perm == 0][model_features == c("PhyloRan_LocRan")],
    aes(x=chr, color =chr ,y=(N)), size = 3) +
    facet_wrap(~bin.id, scales = "free_y") ->
    cum.dis.plot
  ggsave(cum.dis.plot, file = "cum.dis.plot.pdf")
