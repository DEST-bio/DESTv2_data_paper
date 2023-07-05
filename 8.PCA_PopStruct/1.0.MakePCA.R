###
###

### libraries
library(FactoMineR)
library(factoextra)
library(SeqArray)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
library(vroom)
library(poolfstat)

####
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_26April2023.csv"
samps <- fread(meta_git)
setDT(samps)

### open GDS
genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds", allow.duplicate=T)

###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim




###
dat.o <- get(load("/project/berglandlab/DEST2.0_working_data/Filtered_30miss/AFmatrix.flt.Rdata"))
samps <- get(load("/project/berglandlab/DEST2.0_working_data/joint.metadata.Rdata"))
#####

grep( "SIM" , samps$sampleId) -> sim.pos
grep( "CN_Bei_Bei_1_1992-09-16" , samps$sampleId) -> Beijing.pos

samps[-c(sim.pos, Beijing.pos),] -> samps.o.nsim.noBeij
#grep( "SIM" , dat.o.nsim.noBeij$sampleId)

### --- filter data ---

dat.o[samps.o.nsim.noBeij$sampleId, 
      #sample(dim(dat.o)[2], 10000) 
      ] -> dat.o.nsim.noBeij

###### Filter SNPs by thinning phisical linakge
source("/home/yey2sn/software/ThinLDinR_SNPtable.R")
# pickSNPs<-function(map, dist = 100000) ...
### --> Working on this <<<<< ----
colnames(dat.o) -> snps.ids

data.frame(SNP_id = snps.ids) %>%
  separate(SNP_id, into = c("chr","pos","id")) ->
  SNPs.df

SNPs.df$pos = as.numeric(SNPs.df$pos)

SNPs.df %>%
  arrange(chr, pos) ->
  SNPs.df.sort

SNPs.df.sort %>% head

### apply thinning
pickSNPs(SNPs.df.sort, dist = 500 ) -> thinned.set
thinned.set %>% length()

SNPs.df.sort[thinned.set,] %>% 
  mutate(SNP_id = paste(chr, pos, id, sep = "_")) ->
  SNPs.df.sort.thinned
  
#########
#########
#########

which(colnames(dat.o.nsim.noBeij) %in%  SNPs.df.sort.thinned$SNP_id) -> SNP.selector

dat.o.nsim.noBeij[,SNP.selector] %>% 
  PCA(graph = F, ncp = 100) ->
  pca.object

save(pca.object, file = "pca.object.Rdata")

#### Plot

pca.object$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) %>%
  filter(!is.na(continent))->
  pca.meta.dim

pca.meta.dim %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = continent
  )) +
  geom_point() +
  theme_bw() ->
  PCA12.dim

#ggsave(PCA12.dim, file = "PCA12.dim.pdf", w = 5, h =4)

pca.meta.dim %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.3,
    color = continent
  )) +
  geom_point() +
  theme_bw() ->
  PCA13.dim

#ggsave(PCA13.dim, file = "PCA13.dim.pdf", w = 5, h =4)
