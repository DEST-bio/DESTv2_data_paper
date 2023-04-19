### Explore Model of DEST 2.0
### 

library(tidyverse)
library(vroom)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(zoo)
library(adegenet)
library(reshape2)
library("geosphere")
library(gmodels)
library(ggforce)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(forcats)

####
####
####
#sftp://rivanna.hpc.virginia.edu/scratch/yey2sn/DEST2_analysis/dapc_dims/xval_province.DEST2.0.Rdata

GIMS = get(load("DEST.2.0.GIMS.Rdata"))
unique(GIMS$SNP_id) -> GIMS.snps.ids

model2.0 <- get(load("xval_province.DEST2.0.Rdata"))
model2.0$DAPC$ind.coord %>% rownames %>% .[-grep("^NA", .)] -> samps.for.model.train
#model2.0$DAPC$pca.loadings %>% rownames() %>% sort -> GIMs
###
message("now loading AF")
AF.d <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
samps %>% filter(sampleId %in% samps.for.model.train) ->
samps.used.in.modeltrain 
####
#### ---> Creating global object 
AF.d[samps.for.model.train,GIMS.snps.ids] -> AF.gims
AF.gims_naImp = na.aggregate(AF.gims)

####
#### Part 1. Split sample set..
#### 
#i=1
out2 = list()
out =
foreach(i=1:dim(AF.gims_naImp)[1],
        #i=1:10,
        .combine = "rbind",
        .errorhandling = "remove"
        )%do%{
          
          message(i)
          
AF.gims_naImp[-which(rownames(AF.gims_naImp) == samps$sampleId[i]),] ->
  AF.gims_naImp_loo
samps.used.in.modeltrain_loo =
  samps.used.in.modeltrain %>% filter(!sampleId == samps$sampleId[i])

AF.gims_naImp[samps$sampleId[i],] ->
  AF.gims_naImp_anchor

##########
dapc(AF.gims_naImp_loo, 
     grp=samps.used.in.modeltrain_loo$province, 
     n.pca=100, n.da=75) ->
  model.lo_out

predict.dapc(model.lo_out, newdata=AF.gims_naImp_anchor) -> LOO_predictions

  data.frame(pop = rownames(t(LOO_predictions$posterior)),
         P = t(LOO_predictions$posterior)) ->
    tmp.p
names(tmp.p)[2] = "P"
  
tmp.p %>%
  slice_max(P) ->
  post.loo

predicted.pop = post.loo$pop
predicted.post = post.loo$P

samps[sampleId == samps$sampleId[i]] %>%
  group_by(province) %>%
  summarize(real.lat = mean(lat),
            real.long = mean(long),
            cont = unique(continent)
            ) -> real.mean

samps[province == post.loo$pop] %>%
  group_by(province) %>%
  summarize(pred.lat = mean(lat),
            pred.long = mean(long),
  ) -> pred.mean

cbind(predicted.pop, 
      predicted.post, 
      pred.mean$pred.lat, 
      pred.mean$pred.long, 
      real.mean, sampleid = samps$sampleId[i]) -> o
out2[[i]] = o
return(o)

        }

save(out, file = "DEST2.0.model.predictions.Rdata")

#### plot lat and long
out %>%
  ggplot(
    aes(
      x=`pred.mean$pred.lat`,
      y=real.lat,
      color = predicted.post)
  ) +
  geom_abline(slope = 1) +
  geom_point() ->
  lat.pred

ggsave(lat.pred, file = "lat.pred.pdf")

out %>%
  ggplot(
    aes(
      x=`pred.mean$pred.long`,
      y=real.long,
      color = predicted.post)
  ) +
  geom_abline(slope = 1) +
  geom_point() ->
  long.pred

ggsave(long.pred, file = "long.pred.pdf")

#### Calculate the haversine distance 


out %>% 
  group_by(sampleid) %>%
  mutate(hav_d = distHaversine(
    matrix(c(`pred.mean$pred.long`, real.long, `pred.mean$pred.lat`, real.lat),nrow = 2))/1000) ->
  out.hvd
  
out.hvd %>%
  ggplot(
    aes(
      x=cont,
      y=(hav_d)
    )
  ) +
  geom_hline(yintercept = 400) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lr")  +
  geom_boxplot() ->
  hav_dist

ggsave(hav_dist, file = "hav_dist.pdf")

out.hvd %>%
  filter(hav_d > 500) %>%
  ggplot(
    aes(
      x=fct_reorder(sampleid, hav_d),
      y=hav_d,
      label=paste(province,predicted.pop,sep = "-->"),
    )
  ) + 
  geom_hline(yintercept = 1000) +
  geom_point() +
  ylim(0, 20000) +
  geom_text(nudge_y = 3000,size = 2) +
  coord_flip() ->
  troubled.samps
  
ggsave(troubled.samps, file = "troubled.samps.pdf")
ggsave(troubled.samps, file = "troubled.samps.png")


