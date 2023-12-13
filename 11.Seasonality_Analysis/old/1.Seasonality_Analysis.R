### Seasonality Analysis


# --> DEBUG SETTINGS
#k=1
#mAF.filt = 0.05
#perms.ith = 0

## model
model = "season_pop_by_year_All"

args = commandArgs(trailingOnly=TRUE)
##########
k = as.numeric(args[1])
##########
perms.ith = as.numeric(args[2])
#
##########
### K --> 1-7966

########
## jobId=1; nJobs=5000

### libraries
library(bigmemory)
library(fastglm)
library(data.table)
library(tidyverse)
#library(gdata)
library(lubridate)
library(foreach)
library(SeqArray)
#library(glmnet)
library(doMC)
registerDoMC(5)
library(magrittr)
library(reshape2)
library(permute)

### Funcs:
anovaFun <- function(m1, m2) {
  #  m1 <- t0; m2<- t1.year
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  
  chisq <- -2*(ll1-ll2)
  
  1-pchisq(chisq, parameter)
  
}

###
###

message("now loading AF")
AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
message("now loading DP")
DP.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/DPmatrix.flt.Rdata"))

### open metadata
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))

### open filtering object
keep.remove <- get(load("/project/berglandlab/DEST2.0_working_data/keep.fail.samps.Rdata"))

pass.filt <- filter(keep.remove, keep == "PASS")$sampleId

### open seasonal pairs
seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))

#seasonal.sets %>% 
#  left_join(samps) %>%
#  group_by(country) %>%
#  summarize(N=n())

set.samps <- filter(seasonal.sets, sampleId %in%  pass.filt)$sampleId

### ---> loop
#SNPs = "2L_10038126_snp1"
#AF.d %>% colnames %>% head -> tester

message("Preparing iterations")

#PREPARE SNP ITERATORs
dim(AF.d)[2] -> total.snps
snp.indexes = 1:total.snps
colnames(AF.d) -> SNPguides

guide = data.table(index = snp.indexes, SNP = SNPguides)

### define windows
win.snp <- 500
step.snp <- 501

## prepare windows
MASTER_GUIDE = 
data.table(start=seq(from=min(guide$index), 
                     to=max(guide$index)-win.snp, by=step.snp),
           end=seq(from=min(guide$index), 
                   to=max(guide$index)-win.snp, 
                   by=step.snp) + win.snp)
           
SELECTED_SNPSs = SNPguides[MASTER_GUIDE$start[k]:MASTER_GUIDE$end[k]]

message("Launching Iterations")
#####
#####
#####
#####
#####
o.glm =
foreach(SNPs = SELECTED_SNPSs,
        .combine = "rbind",
        .errorhandling = "remove"
        )%do%{

            
AF.tmp = AF.d[set.samps, SNPs]
DP.tmp = DP.d[set.samps, SNPs]

##
# How much missing data?
Nsamps = length(AF.tmp)
Nmiss = sum(is.na(AF.tmp))
N0 = sum(AF.tmp == 0, na.rm = T)
N1 = sum(AF.tmp == 1, na.rm = T)
Mean_AF = mean(AF.tmp, na.rm = T)
Median_AF = median(AF.tmp, na.rm = T)
Mean_DP = mean(DP.tmp, na.rm = T)

###
# --> Make DF 
data.frame(
  sampleId = set.samps,
  AF = AF.tmp,
  DP = DP.tmp
) %>% left_join(dplyr::select(samps, sampleId, nFlies, locality, year, set )) %>% 
  mutate(nEff=round((DP*2*nFlies - 1)/(DP+2*nFlies))) %>%
  mutate(af_nEff=round(AF*nEff)/nEff) %>%
  left_join(dplyr::select(seasonal.sets, season))->
  tmp.af

Mean_NEFF = mean(tmp.af$nEff, na.rm = T)


tmp.af %>%
  filter(!is.na(AF)) %>%
  group_by(locality,year,season) %>%
  slice_head() %>%
  reshape2::dcast(locality+year ~ season, value.var = "sampleId") %>%
  .[complete.cases(.),] -> approved_locs.df

approved_locs.vec = c(approved_locs.df$fall, approved_locs.df$spring)
tmp.af %<>% filter(sampleId %in% approved_locs.vec)

### Check consistentcy of pairs
dcast(tmp.af, locality+year ~ season, value.var = "af_nEff") %>%
  mutate(deltaAF =  fall-spring) -> deltaAF.df
###
TotalPairs = dim(deltaAF.df)[1]
MeanDelta = mean(abs(deltaAF.df$deltaAF), na.rm = T)
Neg_pairs = sum(deltaAF.df$deltaAF < 0, na.rm = T)
Pos_pairs = sum(deltaAF.df$deltaAF > 0, na.rm = T)
Zero_pairs = sum(deltaAF.df$deltaAF == 0, na.rm = T)
  
###
tmp.af$season = as.factor(tmp.af$season)
tmp.af$locality = as.factor(tmp.af$locality)
tmp.af$year = as.factor(tmp.af$year)
 
#mean(tmp.af$AF) -> mean.AF

#if(mean.AF >= mAF.filt){
  
  glm.method <- 0
  
  y <- tmp.af$af_nEff
  X.loc.year <- model.matrix(~locality:year, tmp.af)
  X.seas.loc.year <- model.matrix(~locality:year+season, tmp.af)
  
  t.H0 <- fastglm(x=X.loc.year, 
                  y=y, 
                  family=binomial(), 
                  weights=tmp.af$nEff, method=glm.method)
  
  t.H1 <- fastglm(x=X.seas.loc.year, 
                  y=y, 
                  family=binomial(), 
                  weights=tmp.af$nEff, method=glm.method)
  
  ####
  chr = strsplit(SNPs, "_")[[1]][1]
  pos = strsplit(SNPs, "_")[[1]][2]
  
  data.table(
    chr = chr,
    pos = pos,
    perm = 0,
    AIC=c(AIC(t.H0)),
    b_season=t.H1$coef[2], 
    se_season=t.H1$se[2],
    nObs_i=dim(tmp.af)[1],
    p_lrt=anovaFun(t.H0, t.H1),
    ##Extra info
    SNP_id=SNPs,
    Nsamps_tot=Nsamps,
    Nmiss=Nmiss,
    N0=N0,
    N1=N1,
    Mean_AF=Mean_AF,
    Median_AF=Median_AF,
    Mean_DP=Mean_DP,
    Mean_NEFF=Mean_NEFF,
    TotalPairs=TotalPairs,
    MeanDelta=MeanDelta,
    Neg_pairs=Neg_pairs,
    Pos_pairs=Pos_pairs,
    Zero_pairs=Zero_pairs,
    model=model
  ) -> obs.data
  
  obs.data %<>% mutate(N_pairs_loss = nObs_i/Nsamps )
  
  #### do permutations
  if(perms.ith > 0){
    perms.df =
      foreach(perm = 1:perms.ith, 
              .combine = "rbind", 
              .errorhandling = "remove")%do%{
                
                set.seed(perm)
                
                glm.method <- 0
                
                suffle_vector <- shuffle(tmp.af$af_nEff)
                
                y.shuff <- tmp.af$af_nEff[suffle_vector]
                X.loc.year <- model.matrix(~locality:year, tmp.af)
                X.seas.loc.year <- model.matrix(~locality:year+season, tmp.af)
                
                #tmp.af$nEff[suffle_vector]*y.shuff
                
                t.H0.p <- fastglm(x=X.loc.year, 
                                  y=y.shuff, 
                                  family=binomial(), 
                                  weights=tmp.af$nEff[suffle_vector], method=glm.method)
                
                t.H1.p <- fastglm(x=X.seas.loc.year, 
                                  y=y.shuff, 
                                  family=binomial(), 
                                  weights=tmp.af$nEff[suffle_vector], method=glm.method)
                
                ####
                data.table(
                  chr = chr,
                  pos = pos,
                  perm = perm,
                  AIC=c(AIC(t.H0.p)),
                  b_season=t.H1.p$coef[2], 
                  se_season=t.H1.p$se[2],
                  nObs_i=dim(tmp.af)[1],
                  p_lrt=anovaFun(t.H0.p, t.H1.p),
                  ##Extra info
                  SNP_id=SNPs,
                  Nsamps_tot=Nsamps,
                  Nmiss=Nmiss,
                  N0=N0,
                  N1=N1,
                  Mean_AF=Mean_AF,
                  Median_AF=Median_AF,
                  Mean_DP=Mean_DP,
                  Mean_NEFF=Mean_NEFF,
                  TotalPairs=TotalPairs,
                  MeanDelta=MeanDelta,
                  Neg_pairs=Neg_pairs,
                  Pos_pairs=Pos_pairs,
                  Zero_pairs=Zero_pairs,
                  model=model
                )  -> o.inner
                o.inner %<>% mutate(N_pairs_loss = nObs_i/Nsamps )
                return(o.inner)
              } ## inner
  } ### PERMS
  if(perms.ith <= 0){
    message("No Permutatations done")
    perms=c()}
    
return(rbind(obs.data, perms.df))

} ## outer loop

#####
#####
#####

#o.glm %>% 
#  filter(nObs_i / nObs_tot > 0.8 ) %>%
#  ggplot(aes(
#    p_lrt
#  )) + 
#  geom_histogram() ->
#  hist.test
#ggsave(hist.test, file = "hist.test.pdf")

#####
#####

name.of.file = paste("GLM",model,
                     "SNPrange",MASTER_GUIDE$start[k],MASTER_GUIDE$end[k], "Rdata", 
                     sep ="." )

path.to.save = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_out.03.31.2023"
full.path = paste(path.to.save,name.of.file, sep = "/")

save(o.glm,
     file = full.path)

