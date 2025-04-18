library(foreach)
library(data.table)
library(gmodels)
library(foreach)
library(tidyverse)
library(magrittr)
library(broom)

###
fst_info <- fread("/gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/FST.14900.txt")
mean(fst_info$V1)
mean(fst_info$V2)
mean(fst_info$V3)
###
files <- system(paste("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/ag.results"),
                intern = T)

o.ag =
  foreach(i = files,
          .combine = "rbind")%do%{
    
            tmp <- get(load(
              paste("/gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/ag.results/",i, sep = "")
              ))
          }

o.ag %>%
  group_by(replicate) %>%
  summarise(minE = min(Estimate),
            maxE = max(Estimate)) %>%
  filter(minE > 0 & maxE < 1) %>%
  .$replicate -> 
  pass.filt

o.ag %<>% filter(replicate %in% pass.filt)

###
o.ag %>%
  group_by(replicate) %>%
  summarize(cor=cor(Estimate, ancestry)) %>%
  ggplot(aes(
    cor)) +
  geom_histogram() +
  theme_bw() ->
  corr.plot

ggsave(corr.plot, file = "corr.plot.pdf",
       w=4, h=4)

cor.test(o.ag$Estimate, o.ag$ancestry)
###

o.ag %>%
  dplyr::select(Estimate, ancestry, sampleId, replicate) %>%
  melt(id = c("sampleId", "replicate")) %>%
  ggplot(aes(
    x=as.character(sampleId),
    y=value,
    fill=variable)) +
  geom_boxplot() +
  geom_hline(yintercept =  0, linetype = "dashed") +
  geom_hline(yintercept =  1, linetype = "dashed") +
  theme_bw() ->
  estimates.plot

ggsave(estimates.plot, file = "estimates.plot.pdf",
       w=6,h=4)

o.ag %>%
  mutate(diff = (Estimate-ancestry)) %>%
  ggplot(aes(
    x=as.character(sampleId),
    y=diff)) +
  geom_hline(yintercept =  0.5, linetype = "dashed") +
  geom_hline(yintercept =  -0.5, linetype = "dashed") +
  geom_hline(yintercept =  0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() ->
  diff.plot

ggsave(diff.plot, file = "diff.plot.pdf",
       h=4, w=4)


o.ag %>%
  mutate(diff = (Estimate-ancestry)) %>%
  group_by(sampleId) %>%
  summarise(med = median(diff))

###

o.ag %>%
  dplyr::select(Estimate, ancestry, sampleId, replicate) %>%
  melt(id = c("sampleId", "replicate")) %>%
  ggplot(aes(
    x=(sampleId),
    y=value,
    group=replicate,
    color=variable)) +
  facet_grid(~variable) +
  geom_smooth(method = "lm", se = F, alpha = 0.3) +
  theme_bw()  ->
  traject.plot

ggsave(traject.plot, file = "traject.plot.pdf",
       w=6, h=4)

##
o.ag %>%
  nest_by(replicate) %>%
  mutate(mod = list(lm(Estimate ~ sampleId, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "sampleId") -> EST_beta
names(EST_beta)[c(3,4)] = c("estimate.ES" , "std.error.ES")

o.ag %>%
  nest_by(replicate) %>%
  mutate(mod = list(lm(ancestry ~ sampleId, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "sampleId") -> ANS_beta
names(ANS_beta)[c(3,4)] = c("estimate.AN" , "std.error.AN")

cor.test(ANS_beta$estimate.AN, EST_beta$estimate.ES)

left_join(
ANS_beta[,c("replicate" ,"estimate.AN" , "std.error.AN")],
EST_beta[,c("replicate" ,"estimate.ES" , "std.error.ES")]) %>%
  ggplot(aes(
    x=estimate.AN,
    y=estimate.ES,
    xmin=estimate.AN-std.error.AN,
    xmax=estimate.AN+std.error.AN,
    ymin=estimate.ES-std.error.ES,
    ymax=estimate.ES+std.error.ES,
  )) + geom_point(alpha = 0.2) + 
  geom_errorbar(alpha = 0.2) + 
  geom_errorbarh(alpha = 0.2) + 
  geom_density2d() + theme_bw() +
  geom_abline(slope = 1)->
  betas_plot

ggsave(betas_plot, file = "beta.plot.pdf",
       w=4, h=4)

left_join(
  ANS_beta[,c("replicate" ,"estimate.AN" , "std.error.AN")],
  EST_beta[,c("replicate" ,"estimate.ES" , "std.error.ES")]) %>%
  dplyr::select(replicate, estimate.AN, estimate.ES) %>%
  melt(id = "replicate") %>%
  ggplot(aes(
    x=variable,
    y=value
  )) + geom_boxplot() + 
  theme_bw() ->
  betas_boxplot

ggsave(betas_boxplot, file = "betas_boxplot.pdf",
       h=4,w=4)


