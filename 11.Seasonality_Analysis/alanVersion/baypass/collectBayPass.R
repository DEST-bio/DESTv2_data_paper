# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(dplyr)

### SNPs
   snps <- fread("/standard/vol186/bergland-lab/Gio/subbaypass/dest_pos_table.txt")
   setnames(snps, "subpoolMRK", "MRK")
   setkey(snps, chr, pos)

### load SNP definition files
  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/", "snpdet", full.name=T)
  snp.dt <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    file <- tstrsplit(fl.i, "/")%>%last()
    subPool <- gsub(".snpdet", "", tstrsplit(file, "_")[[2]])%>%as.numeric()

    tmp <- fread(fl.i)
    tmp[,subpool:=subPool]
    tmp[,MRK:=1:dim(tmp)[1]]
    tmp
  }
  snp.dt <- rbindlist(snp.dt)
  setnames(snp.dt, c("V1", "V2", "V3", "V4"), c("chr", "pos", "ref", "alt"))
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snp.dt, snps[,-c("MRK", "subpool"),with=F])

### get files
  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass", "xtx", full.name=T)

### iterate
  xtx <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    file <- tstrsplit(fl.i, "/")%>%last()
    subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
    rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()

    tmp <- fread(fl.i)
    setnames(tmp, "log10(1/pval)", "negLog10P")
    tmp[,subpool:=subPool]
    tmp[,rep:=rep]
    tmp
  }
  xtx <- rbindlist(xtx)

## merge
  setkey(snp.dt, MRK, subpool)
  setkey(xtx, MRK, subpool)
  xtx <- merge(xtx, snp.dt)
  save(xtx, file="/scratch/aob2x/xtx_dest2.raw.Rdata")

### average
  xtx.ag <- xtx[,list(XtXst_median=median(XtXst), XtXst_mean=mean(XtXst),
                      neglogp_median=median(negLog10P), neglogp_mean=mean(negLog10P)), list(MRK, subpool, chr, pos, invName)]
  xtx.ag[order(subpool)]
  table(xtx.ag$XtXst_mean>200, xtx.ag$chr, xtx.ag$invName!="noInv")

##########
### contrast
##########

### get files
  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/", "contrast", full.name=T)

  cont <- foreach(fl.i=fl)%dopar%{
    message(fl.i)
    # fl.i <- fl[1]
    file <- tstrsplit(fl.i, "/")%>%last()
    subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
    rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()

    tmp <- fread(fl.i)
    setnames(tmp, "log10(1/pval)", "negLog10P")
    tmp[,subpool:=subPool]
    tmp[,rep:=rep]
    tmp
  }
  cont <- rbindlist(cont)

## merge
  setkey(snp.dt, MRK, subpool)
  setkey(cont, MRK, subpool)
  cont <- merge(cont, snp.dt)
  save(cont, file="/scratch/aob2x/cont_dest2.raw.Rdata")

### average
  cont.ag <- cont[,list(C2_std_median=median(C2_std), C2_std_mean=mean(C2_std),
                        neglogp_median=median(negLog10P), neglogp_mean=mean(negLog10P)), list(MRK, subpool, chr, pos)]

  cont.ag[,p:=10^(-neglogp_mean)]
  cont.ag[,q:=p.adjust(p, "fdr")]
  table(cont.ag$q<.15, cont.ag$chr)

### GLM
  load(file="/standard/vol186/bergland-lab/DEST_v2/GLMER_output_EuropeSeasonality.Rdata")
  o[,rnp:=rank(p_lrt)/(1+length(p_lrt))]
  setkey(o, chr, pos)
  setkey(xtx.ag, chr, pos)
  setkey(cont.ag, chr, pos)

  m <- merge(o, xtx.ag)
  m <- merge(m, cont.ag)
  m[,q:=p.adjust(p_lrt, "fdr")]

### save
  save(m, file="~/dest2_glm_baypass.Rdata")
  fisher.test(table(m$rnp<.05, m$XtXst_median>450))
  fisher.test(table(m$rnp<.05, m$C2_std_median>5))
  fisher.test(table(m$q.x<.0005, m$q.y<.05))

  fisher.test(table(m$q.x<.0005 & m$q.y<.05, m$q<.05))

  m[which.max(C2_std)]


### local
  library(data.table)
  library(ggplot2)
  library(qvalue)

### data
  load("~/dest2_glm_baypass.Rdata")
  m[,cont_p:=p]
  m[,cont_q:=qvalue(cont_p)]

### basic plot
  ggplot(data=m, aes(x=-log10(p_lrt), y=)) + geom_hex() + facet_grid(I(invName.x!="noInv") ~ chr)
  ggplot(data=m, aes(x=-log10(p_lrt), y=C2_std_median)) + geom_hex() + facet_grid(I(invName.x!="noInv") ~ chr)
  ggplot(data=m, aes(x=XtXst_median, y=C2_std_median)) + geom_hex() + facet_grid(I(invName.x!="noInv") ~ chr)

  ggplot(data=m, aes(XtXst_mean)) + geom_histogram() + facet_grid(I(invName.x!="noInv") ~ chr)

  ggplot(data=m[XtXst_median>200], aes(y=XtXst_median, x=pos, color=invName.x!="noInv")) +
  geom_point() + facet_grid( ~ chr)

  ggplot(data=m[C2_std_median>5], aes(y=C2_std_mean, x=pos, color=invName.x)) +
  geom_point(size=.5, alpha=.5) +
  facet_grid( ~ chr)

m[C2_std_median>20][chr=="3R"]

  buffer <- 10000
  table(m$chr=="3R" & m$pos>=(13215951-buffer) & m$pos<=(13269700+buffer),
        m$C2_std_median>15)%>%fisher.test()


  fisher.test(table(m[chr=="2R"]$XtXst_median>350, m[chr=="2R"]$rnp<.05))
