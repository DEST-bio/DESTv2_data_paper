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

### Bring - MOMENTS
filesM <- system(paste("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/moments_slim_output"),
                intern = T)

o.MO =
  foreach(i = filesM,
          .combine = "rbind", .errorhandling = "remove")%do%{
            
            message(i)
            
            tmp <- fread(
              paste("/gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/moments_slim_output/",i, sep = "")
            ) 
            
            tmp %>%
              group_by(V1) %>%
              summarise(MO.Admx = mean(admix_prop)) %>%
              separate(V1, into = c("etc1","sampleId", "etc2","replicate")) %>%
              dplyr::select(sampleId, replicate, MO.Admx)
            
          }

o.ag$sampleId = as.numeric(o.ag$sampleId)
o.MO$sampleId = as.numeric(o.MO$sampleId)

left_join(o.ag, o.MO) ->
  joint.o.MO.LT

save(joint.o.MO.LT, file = "simulated.SLIM.moments.LM.Rdata")
load("simulated.SLIM.moments.LM.Rdata")

###
joint.o.MO.LT %>%
  dplyr::select(sampleId, ancestry, replicate , Estimate,  MO.Admx) %>%
  melt(id = c("sampleId", "ancestry", "replicate" )) %>%
  group_by(replicate, variable) %>%
  summarize(cor=cor(value, ancestry)) %>%
  ggplot(aes(
    cor, fill = variable)) +
  geom_histogram() +
  theme_bw() ->
  corr.plot

ggsave(corr.plot, file = "corr.plot.pdf",
       w=4, h=4)

cor.test(joint.o.MO.LT$Estimate, joint.o.MO.LT$ancestry)
#0.7975302 p-value < 2.2e-16
cor.test(joint.o.MO.LT$MO.Admx, joint.o.MO.LT$ancestry)
#0.8013978p-value < 2.2e-16

joint.o.MO.LT %>%
  dplyr::select(sampleId, ancestry, replicate , Estimate,  MO.Admx) %>%
  melt(id = c("sampleId", "replicate" )) %>% 
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

joint.o.MO.LT %>%
  mutate(diff_LM = (Estimate-ancestry),
         diff_MO = (MO.Admx-ancestry)) %>%
  dplyr::select(sampleId,diff_LM,diff_MO ) %>%
  melt(id = "sampleId") %>%
  ggplot(aes(
    x=as.character(sampleId),
    y=value,
    fill = variable)) +
  geom_hline(yintercept =  0.5, linetype = "dashed") +
  geom_hline(yintercept =  -0.5, linetype = "dashed") +
  geom_hline(yintercept =  0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() ->
  diff.plot

ggsave(diff.plot, file = "diff.plot.pdf",
       h=4, w=4)


joint.o.MO.LT %>%
  mutate(diff_LM = as.numeric(Estimate-ancestry),
         diff_MO = as.numeric(MO.Admx-ancestry)) %>%
  dplyr::select(sampleId,diff_LM,diff_MO ) %>%
  melt(id = "sampleId") %>%
  group_by(variable, sampleId) %>% 
  summarise(m.var = mean(value, na.rm = T))

###

joint.o.MO.LT %>%
  dplyr::select(sampleId, ancestry, replicate , Estimate,  MO.Admx) %>%
  melt(id = c("sampleId", "replicate" )) %>% 
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
joint.o.MO.LT %>%
  nest_by(replicate) %>%
  mutate(mod = list(lm(Estimate ~ sampleId, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "sampleId") -> EST_beta
names(EST_beta)[c(3,4)] = c("estimate.ES" , "std.error.ES")

joint.o.MO.LT %>%
  filter(!is.na(MO.Admx)) %>%
  nest_by(replicate) %>%
  mutate(mod = list(lm(MO.Admx ~ sampleId, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "sampleId") -> MO_beta
names(MO_beta)[c(3,4)] = c("estimate.MO" , "std.error.MO")

joint.o.MO.LT %>%
  nest_by(replicate) %>%
  mutate(mod = list(lm(ancestry ~ sampleId, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "sampleId") -> ANS_beta
names(ANS_beta)[c(3,4)] = c("estimate.AN" , "std.error.AN")


left_join(
ANS_beta[,c("replicate" ,"estimate.AN" , "std.error.AN")],
EST_beta[,c("replicate" ,"estimate.ES" , "std.error.ES")]) %>%
left_join(MO_beta[,c("replicate" ,"estimate.MO" , "std.error.MO")]) ->
  joint.slopes

cor.test(joint.slopes$estimate.AN, joint.slopes$estimate.ES)
cor.test(joint.slopes$estimate.AN, joint.slopes$estimate.MO)

joint.slopes %>%
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

joint.slopes %>%
  ggplot(aes(
    x=estimate.AN,
    y=estimate.MO,
    xmin=estimate.AN-std.error.AN,
    xmax=estimate.AN+std.error.AN,
    ymin=estimate.ES-std.error.MO,
    ymax=estimate.ES+std.error.MO,
  )) + geom_point(alpha = 0.2) + 
  geom_errorbar(alpha = 0.2) + 
  geom_errorbarh(alpha = 0.2) + 
  geom_density2d() + theme_bw() +
  geom_abline(slope = 1)->
  betas_plot.MO

ggsave(betas_plot.MO, file = "betas_plot.MO.pdf",
       w=4, h=4)



joint.slopes %>%
  dplyr::select(replicate, estimate.AN, estimate.ES, estimate.MO) %>%
  melt(id = "replicate") %>%
  ggplot(aes(
    x=variable,
    y=value
  )) + geom_boxplot() + 
  theme_bw() ->
  betas_boxplot

ggsave(betas_boxplot, file = "betas_boxplot.pdf",
       h=4,w=4)

####

