
require(poolfstat)
library(SeqArray)
library(foreach)
library(poolfstat)
library(data.table)
library(tidyverse)
library(magrittr)

i=34


genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")
meta <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")
guide <- fread("/gpfs2/scratch/jcnunez/DEST2.0_analysis/MOMENTS_REVISION/pairs_guide_file.txt")

message("making snp table")
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")
## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

  
samps <- c(guide$Parent1[i],
           guide$Parent2[i],
           guide$Derived[i],
           "SIM_SIM_w501_1_NA-MM-DD",
           "US_Cal_Sto_1_2013-09-03"
           )

samps <-  sort(samps)

seqResetFilter(genofile)
seqSetFilter(genofile, 
             sample.id=samps,
             variant.id=filter(snps.dt, chr %in%
                                 c("2L", "2R", "3L", "3R"))$variant.id)

filter(meta, sampleId %in% samps)$nFlies
#seqGetData(genofile, "sample.id")

ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

cbind(snps.dt, t(dp), t(ad$data)) %>%
  .[complete.cases(.),] ->
  dpt

tmp.sim=which(seqGetData(genofile, "sample.id")=="SIM_SIM_w501_1_NA-MM-DD")

MinDP <- apply(dpt[,c(8:12)[-tmp.sim]], 1, function(x) min(x))
FixL <- apply(dpt[,c(13:17)[-tmp.sim]], 1, function(x) sum(x == 0))


dpt[c("chr", "pos", "variant.id", "ref_allele", "alt_allele")] %>% mutate(MinDP = MinDP,
                FixL = FixL) %>% 
  filter(MinDP > 5) %>%
  filter(FixL < 3) -> dpt

seqSetFilter(genofile, 
             sample.id=samps,
             variant.id=dpt$variant.id)

ad2 <- seqGetData(genofile, "annotation/format/AD")
dp2 <- seqGetData(genofile, "annotation/format/DP")


pool <- new("pooldata",
            npools=dim(ad2$data)[1], #### Rows = Number of pools
            nsnp=dim(ad2$data)[2], ### Columns = Number of SNPs
            refallele.readcount=t(ad2$data),
            readcoverage=t(dp2),
            poolsizes=filter(meta, sampleId %in% samps)$nFlies,
            poolnames = samps,
            snp.info = dpt[,c("chr","pos","ref_allele","alt_allele")])


gc()
tmp.sim=which(pool@poolnames=="SIM_SIM_w501_1_NA-MM-DD")
pool@readcoverage[,tmp.sim]=pool@readcoverage[,tmp.sim]*100
pool@refallele.readcount[,tmp.sim]=pool@refallele.readcount[,tmp.sim]*100
gc()

tmp.sel=pooldata.subset(pool,
                        pool.index = 
                          which(pool@poolnames%in%
                                  c(guide$Parent1[i],"SIM_SIM_w501_1_NA-MM-DD",guide$Parent2[i],guide$Derived[i],"US_Cal_Sto_1_2013-09-03")),
                        min.maf = 0.05, min.cov.per.pool = 5)

tmp.sel.f4=compute.fstats(pool,
                          nsnp.per.bjack.block = 10000,
                          computeF4 = T, 
                          return.F2.blockjackknife.samples = TRUE)

f4.ratio=compute.f4ratio(tmp.sel.f4,
                         num.quadruplet = paste0(guide$Parent1[i],",SIM_SIM_w501_1_NA-MM-DD;",guide$Parent1[i],",",guide$Derived[i]),
                         den.quadruplet = paste0(guide$Parent1[i],",SIM_SIM_w501_1_NA-MM-DD;",guide$Parent1[i],",US_Cal_Sto_1_2013-09-03")
                         )
if(i==1){res.f4.ratio=f4.ratio}else{res.f4.ratio=rbind(res.f4.ratio,f4.ratio)}



if(i==1){res.f4.ratio=f4.ratio}else{res.f4.ratio=rbind(res.f4.ratio,f4.ratio)}
cat(tmp.sel.ena[i],"\n")
}

rownames(res.f4.ratio)=tmp.sel.ena
res.f4.ratio=cbind(res.f4.ratio,dmel.info[tmp.sel.ena,]$lat)
res.f4.ratio[,1:5]=1-res.f4.ratio[,1:5] #to estimate African ancestry
write.table(file="res.f4.ratio",res.f4.ratio,quote=F)

pdf("F4ratioEstimates.pdf",h=9,w=9)
plot(res.f4.ratio[,6],res.f4.ratio[,1],xlab="latitude",ylab="Prop. Of African Ancestry",pch=16,
     main="F4-ratio estimates of African ancestry in North-Eastern American samples")
text(res.f4.ratio[,6],res.f4.ratio[,1],tmp.sel.ena,cex=0.5)
graphics.off()
