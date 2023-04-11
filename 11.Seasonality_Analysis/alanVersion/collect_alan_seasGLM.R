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
files.seas = system("ls ./GLM_ALAN_APR102023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_ALAN_APR102023/"

#### binning analysis
seas.p.bin = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/")))
    
    return(tmp)
    
  }

#### Collect into permutation sets
outfolder = "/project/berglandlab/jcbnunez/Shared_w_Alan/GLM_ALAN_APR102023_by_PERM_seas"

### save as independent files
  foreach(p = 1:100)%do%{
    
    message(paste(p, 
                  sep = "/"))
    
    dat.flt %>% filter(perm == p 
    ) -> dat.p
    message(p)
    save(dat.p,
         file = paste(outfolder,
                      paste("GLM_out.perm_id", p, "SeasAlan", "Rdata", sep = "."),
                      sep = "/"
                      ))
  }



#### Analyze the p-value
seas.p.bin ->
dat.flt

ot_across_perms=
foreach(p = 0:100, .combine = "rbind" )%do%{
  #foreach(chr = c("2L","2R","3L","3R"), .combine = "rbind")%do%{
    message(paste(p, 
                  #chr, 
                  sep = "/"))
    
    dat.flt %>% filter(perm == p 
                       #& chr == chr
                       ) -> dat.0
    hist(dat.0$p_lrt, breaks = 80) -> hist.dat
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

  ggplot() +
  geom_violin(
    data = ot_across_perms[perm != 0],
    aes(x= as.factor(hist.dat.mids),
        y=hist.dat.counts),
    fill = "grey", size = 0.5
  ) +
    geom_point(
      data = ot_across_perms[perm == 0],
      aes(x= as.factor(hist.dat.mids),
          y=hist.dat.counts),
      size = 1,
      color = "red"
    ) + 
    #facet_wrap(~chr) +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) -> p.vals.filts
  
  ggsave(p.vals.filts, file = "p.vals.filts.pdf", w = 9, h = 4)
  
    
  
  ggplot(aes(
    x=hist.dat.mids,
    y=(hist.dat.counts),
    size = perm==0,
    color = perm==0,
    alpha = perm==0,
    group = perm
  )) + geom_line() +
  scale_size_manual(values = c(0.5, 2)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~chr, nrow = 1) -> p.vals.filts

ggsave(p.vals.filts, file = "p.vals.filts.pdf", w = 8, h = 2.3)






####
####
####
####

### --> Filter mising 0.9 ; mAF 0.2

perc_dat.ith = seq(from = 0.1, to = 0.9, by = 0.1)
mAF.snp.ith = seq(from = 0.1, to = 0.9, by = 0.1)
seas.p.bin %<>%
  mutate(Zero.snp = N0/Nsamps_tot)

p.vals.range = 
foreach(k = perc_dat.ith, .combine = "rbind")%do%{
  o.2 =
  foreach(i = mAF.snp.ith, .combine = "rbind")%do%{
   
    message(paste(k,i, sep = "|"))
    seas.p.bin %>%
    filter(perm == 0 & N_pairs_loss > k &  Zero.snp > i) -> dat.0
    
    hist(dat.0$p_lrt) -> hist.dat
    data.frame(hist.dat$mids ,    hist.dat$counts) %>%
      mutate(N_pairs_loss = k,
             Zero.snp = i
             ) -> ot
    
    return(ot)} 
  return(o.2)}

p.vals.range %>%
  filter(N_pairs_loss == 0.1) %>% 
  ggplot(aes(
    x=hist.dat.mids,
    y=(hist.dat.counts),
    color = N_pairs_loss,
    group = N_pairs_loss
  )) + geom_line() +
  geom_point() +
  facet_wrap(N_pairs_loss~Zero.snp, nrow = 2, scales = "free_y") -> p.vals.filts

ggsave(p.vals.filts, file = "p.vals.filts.pdf")

###

seas.p.bin %>%
filter(N_pairs_loss < 0.7) %>% 
  filter(perm == 0.0) %>% 
  filter(Zero.snp < 0.7) %>% 
  ggplot(aes(
    x=p_lrt,
    y=b_season
  )) +
  geom_point() ->
  beta_p

ggsave(beta_p, file = "beta_p.png")
  
### 
seas.p.bin %>%
  mutate(beta_tag = case_when(abs(b_season) > 5 ~ "abonrmal",
                              TRUE ~ "normal")) ->
  seas.p.bin

####
seas.p.bin %>%
  filter(N_pairs_loss < 0.7) %>% 
  filter(perm == 0.0) %>% 
  filter(Zero.snp < 0.7) %>% 
  ggplot(aes(
    x=p_lrt,
    y=b_season,
    shape =beta_tag,
    color = N0
  )) +
  geom_point() ->
  beta_p

ggsave(beta_p, file = "beta_p.png")

### N0,, 

seas.p.bin %>%
  dplyr::select(Nmiss, N0, Mean_DP, beta_tag) %>% 
  melt(id = "beta_tag") %>% 
  group_by(variable, beta_tag) %>%
  summarise(m = mean(value))
  
  ggplot(
    aes(
      x=beta_tag,
      y=log10(value)
    )
  ) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") ->
  abnormal.sts

ggsave(abnormal.sts, file = "abnormal.sts.png")




### 
seas.p.bin %>%
  ggplot(aes(
    x=hist.dat.mids,
    y=(hist.dat.counts),
    group = N_pairs_loss,
    color =beta_tag
  )) + geom_line() +
  geom_point() +
  facet_wrap(N_pairs_loss~Zero.snp, nrow = 2, scales = "free_y") -> p.vals.filts

ggsave(p.vals.filts, file = "p.vals.filts.pdf")



###
p.vals.range %>% 
  group_by(perc_dat) %>%
  summarize(Ns = sum(hist.dat.counts))




#######

seas.p.bin %>%
filter(N0 < 50 & beta_tag == "normal") -> dat.flt

hist.p = 
foreach(perm.th = 0:100, .combine = "rbind" )%do%{
  message(perm.th)
  dat.flt %>% filter(perm == perm.th) -> tmp0
    hist(tmp0$p_lrt) -> hist.dat
    data.frame(hist.dat$mids ,    hist.dat$counts) %>%
      mutate(perm = perm.th
      ) -> ot
    return(ot)
  }

hist.p %>%
  ggplot(aes(
    x=hist.dat.mids,
    y=(hist.dat.counts),
    color = perm == 0,
  )) + geom_point()  -> p.vals.filts.o.p

ggsave(p.vals.filts.o.p, file = "p.vals.filts.o.p.pdf")


#######
####  
####  
####  
####  
####  
  seas.p.bin %>%
    filter(perm == 0) %>%
    filter(p_lrt > 0.98) ->
    seas.p.bin.p1
  
  message("now loading AF")
  AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
  
  colnames(AF.d) -> AF.snp.info
  
  data.frame(SNP_id = AF.snp.info) %>%
    separate(SNP_id, into = c("chr", "pos", "snp"),
             remove = F) -> SNPinfo
  
  SNPinfo %>% 
    filter(chr == "2L") ->
    SNPinfo.2L
  
  SNPinfo.2L %>%
    filter(pos %in% seas.p.bin.p1$pos) ->
    SNPinfo.2L.glm
  
####
  AF.d[set.samps,SNPinfo.2L.glm$SNP_id ] %>% 
    t() %>% melt() -> AF.d.samps.sel
  
  names(AF.d.samps.sel)[2] = "sampleId"
  
  AF.d.samps.sel %>%
    left_join(samps) ->
    AF.d.samps.sel.metas

  AF.d.samps.sel.metas %>%
    filter(Var1 == "2L_10040204_snp145") %>%
    ggplot(aes(
      x=jday,
      y=value,
      #group = city,
      color = as.factor(year)
    )) + 
    geom_line() +
    facet_wrap(~province) ->
    af.p1.plot
  ggsave(af.p1.plot, file = "af.p1.plot.png")
  
  
  
###  
  
seas.p.bin %>%
  group_by(AF_bin,perm,chr) %>%
  summarize(Nsps.ag = sum(Nsp)) ->
  seas.p.bin.ag

save(seas.p.bin.ag, file = "seas.p.bin.ag.Rdata")

#### --> collect object GLMS
#### 
seas.p.glm.2L = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/"))) %>%
      mutate(perc_dat = nObs_i/nObs_tot) %>%
      filter(chr == "2L")
    
    return(tmp)
  }

seas.p.glm.2L$pos = as.numeric(seas.p.glm.2L$pos)
seas.p.glm.2L %<>% arrange(chr, pos) 

save(seas.p.glm.2L, file = "seas.p.glm.2L.Rdata")


