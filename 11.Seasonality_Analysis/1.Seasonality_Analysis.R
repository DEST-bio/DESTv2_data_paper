### Seasonality Analysis

args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])
### K --> 1-7966
model = "season_pop_by_year"

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

### open seasonal pairs
seasonal.sets <- get(load("DEST2.seasonals.plusCore20.flip.met.Rdata"))

#seasonal.sets %>% 
#  left_join(samps) %>%
#  group_by(country) %>%
#  summarize(N=n())

set.samps <- seasonal.sets$sampleId

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
###
# --> Make DF 
data.frame(
  sampleId = set.samps,
  AF = AF.tmp,
  DP = DP.tmp
) %>%
  left_join(dplyr::select(samps, sampleId, nFlies, locality, year, set )) %>% 
  mutate(nEff=round((DP*2*nFlies - 1)/(DP+2*nFlies))) %>%
  mutate(af_nEff=round(AF*nEff)/nEff) %>%
  left_join(dplyr::select(seasonal.sets, season))->
  tmp.af

tmp.af %>%
  filter(!is.na(AF)) %>%
  group_by(locality,year,season) %>%
  slice_head() %>%
  reshape2::dcast(locality+year ~ season, value.var = "sampleId") %>%
  .[complete.cases(.),] -> approved_locs.df

approved_locs.vec = c(approved_locs.df$fall, approved_locs.df$spring)

tmp.af %<>% filter(sampleId %in% approved_locs.vec)

tmp.af$season = as.factor(tmp.af$season)
tmp.af$locality = as.factor(tmp.af$locality)
tmp.af$year = as.factor(tmp.af$year)
 
#tmp.af %>% group_by(set) %>% summarise(N=n())
 
####
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
           nObs_tot=dim(seasonal.sets)[1],
           p_lrt=anovaFun(t.H0, t.H1)
           ) -> obs.data

#### do permutations
perms =
foreach(perm = 1:100, 
        .combine = "rbind", 
        .errorhandling = "remove")%do%{
  
          set.seed(perm)
          
          glm.method <- 0
          
          y.shuff <- tmp.af$af_nEff[shuffle(tmp.af$af_nEff)]
          X.loc.year <- model.matrix(~locality:year, tmp.af)
          X.seas.loc.year <- model.matrix(~locality:year+season, tmp.af)
          
          t.H0.p <- fastglm(x=X.loc.year, 
                          y=y.shuff, 
                          family=binomial(), 
                          weights=tmp.af$nEff, method=glm.method)
          
          t.H1.p <- fastglm(x=X.seas.loc.year, 
                          y=y.shuff, 
                          family=binomial(), 
                          weights=tmp.af$nEff, method=glm.method)
          
          ####
          chr = strsplit(SNPs, "_")[[1]][1]
          pos = strsplit(SNPs, "_")[[1]][2]
          
          data.table(
            chr = chr,
            pos = pos,
            perm = perm,
            AIC=c(AIC(t.H0.p)),
            b_season=t.H1.p$coef[2], 
            se_season=t.H1.p$se[2],
            nObs_i=dim(tmp.af)[1],
            nObs_tot=dim(seasonal.sets)[1],
            p_lrt=anovaFun(t.H0.p, t.H1.p)
          ) 
          
} ## inner

rbind(obs.data, perms)

} ## outer

name.of.file = paste("GLM",model,"SNPrange",MASTER_GUIDE$start[k],MASTER_GUIDE$end[k], "Rdata", 
                     sep ="." )
path.to.save = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_out"
full.path = paste(path.to.save,name.of.file, sep = "/")

save(o.glm,
     file = full.path)

