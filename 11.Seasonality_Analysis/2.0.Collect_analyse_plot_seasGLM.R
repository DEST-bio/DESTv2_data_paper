#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
#####
### open metadata
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
### open filtering object
keep.remove <- get(load("/project/berglandlab/DEST2.0_working_data/keep.fail.samps.Rdata"))
pass.filt <- filter(keep.remove, keep == "PASS")$sampleId
### open seasonal pairs
seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
set.samps <- filter(seasonal.sets, sampleId %in%  pass.filt)$sampleId


#### Central files
files.seas = system("ls ./GLM_out.03.27.2023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_out.03.27.2023"

#### binning analysis

seas.p.bin = 
  foreach(fil = files.seas, .combine = "rbind" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/"))) %>%
      mutate(perc_dat = nObs_i/nObs_tot)

    #perform binning with custom breaks
    #tmp %>% 
    #  filter(perc_dat > 0.85) %>%
    #  mutate(AF_bin = chop(p_lrt, breaks=c(
    #  seq(from=0.001, to = 0.009, by = 0.001),
    #  seq(from=0.01, to = 0.09, by = 0.01),
    #  seq(from=0.1, to = 1.0, by = 0.1))
    #  )) -> tmp.i
    #
    #tmp.i %>%
    #  group_by(AF_bin, perm, chr) %>%
    #  summarize(Nsp = n()) %>%
    #  mutate(fil = fil)->
    #tmp.ag
    
    return(tmp)
    
  }

seas.p.bin %>%
  filter(perm == 0) %>%
  ggplot(aes(
    p_lrt
  )) +
  geom_histogram() ->
  hist.test

ggsave(hist.test, file = "hist.test.pdf")

seas.p.bin %>%
  filter(perm == 0 & perc_dat > 0.8) -> dat.0

pdf(file = "beta.hist.test.pdf")
hist(dat.0$p_lrt)
dev.off()

ggsave(beta.hist.test, file = "beta.hist.test.png")

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


