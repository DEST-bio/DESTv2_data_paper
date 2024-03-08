library(data.table)
library(poolfstat)
library(progress)

setwd("/mnt/lustre/scratch/nlsas/home/csic/gcy/mcz/DEST_v2_Fst/ANALYSIS/PCA")

#source("../2compute_hierFST.R")
#source("../generate.jacknife.blocks.R")

# Aut
data_aut=fread(file="../global_fst/aut/out.baypass.geno",data.table=FALSE)
data_2L=fread(file="../global_fst/aut/2L/out.baypass.geno",data.table=FALSE)
data_2R=fread(file="../global_fst/aut/2R/out.baypass.geno",data.table=FALSE)
data_3L=fread(file="../global_fst/aut/3L/out.baypass.geno",data.table=FALSE)
data_3R=fread(file="../global_fst/aut/3R/out.baypass.geno",data.table=FALSE)

data=data_aut
#data=data_2L
#data=data_2R
#data=data_3L
#data=data_3R

data.pool=new("pooldata")
pool.lst=read.csv("../june_metadata.csv")$sampleId
pool.clusters=read.csv("../june_metadata.csv")$cluster2.0_k8
pool.info=read.csv("../june_metadata.csv",h=T,row.names=1)[pool.lst,]
n.pools=length(pool.lst)

data.pool@npools=nrow(pool.info)
idx=seq(1,2*data.pool@npools,2)
data.pool@refallele.readcount=as.matrix(data[,idx])
data.pool@readcoverage=data.pool@refallele.readcount+as.matrix(data[,idx+1])
rm(data) ; gc()
data.pool@nsnp=nrow(data.pool@readcoverage)
data.pool@poolnames=rownames(pool.info)

snp.info_aut=fread("../global_fst/aut/out.snpdet.light.gz",data.table=FALSE)
snp.info_2L=fread("../global_fst/aut/2L/out.snpdet.light.gz",data.table=FALSE)
snp.info_2R=fread("../global_fst/aut/2R/out.snpdet.light.gz",data.table=FALSE)
snp.info_3L=fread("../global_fst/aut/3L/out.snpdet.light.gz",data.table=FALSE)
snp.info_3R=fread("../global_fst/aut/3R/out.snpdet.light.gz",data.table=FALSE)

snp.info=snp.info_aut
#snp.info=snp.info_2L
#snp.info=snp.info_2R
#snp.info=snp.info_3L
#snp.info=snp.info_3R

colnames(snp.info)=c("ID","Chromosome","Position","RefAllele","AltAllele")
rownames(snp.info)=snp.info[,1]
snp.info=snp.info[,-1]
data.pool@snp.info=snp.info
data.pool@poolsizes=2*pool.info$nFlies*2


subset <- pooldata.subset(
  data.pool,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  verbose = TRUE
)

pca_result <- randomallele.pca(subset, scale = TRUE)

pca_result_aut <- pca_result
#pca_result_2L <- pca_result
#pca_result_2R <- pca_result
#pca_result_3L <- pca_result
#pca_result_3R <- pca_result

save(pca_result_aut, file="pca_result_aut.RData")
save(pca_result_2L, file="pca_result_2L.RData")
save(pca_result_2R, file="pca_result_2R.RData")
save(pca_result_3L, file="pca_result_3L.RData")
save(pca_result_3R, file="pca_result_3R.RData")

# X
data=fread(file="../global_fst/X/out.baypass.geno",data.table=FALSE)
data.pool=new("pooldata")
pool.lst=read.csv("../june_metadata.csv")$sampleId
pool.clusters=read.csv("../june_metadata.csv")$cluster2.0_k8
pool.info=read.csv("../june_metadata.csv",h=T,row.names=1)[pool.lst,]

n.pools=length(pool.lst)

data.pool@npools=nrow(pool.info)
idx=seq(1,2*data.pool@npools,2)
data.pool@refallele.readcount=as.matrix(data[,idx])
data.pool@readcoverage=data.pool@refallele.readcount+as.matrix(data[,idx+1])
rm(data) ; gc()
data.pool@nsnp=nrow(data.pool@readcoverage)
data.pool@poolnames=rownames(pool.info)
snp.info=fread("../global_fst/X/out.snpdet.light.gz",data.table=FALSE)
colnames(snp.info)=c("ID","Chromosome","Position","RefAllele","AltAllele")
rownames(snp.info)=snp.info[,1]
snp.info=snp.info[,-1]
data.pool@snp.info=snp.info
data.pool@poolsizes=2*pool.info$nFlies

subset <- pooldata.subset(
  data.pool,
  min.cov.per.pool = -1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  verbose = TRUE
)

pca_result <- randomallele.pca(subset, scale = TRUE)

pca_result_X <- pca_result
save(pca_result_X, file="pca_result_X.RData")
