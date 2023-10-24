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


### get files
  fl <- list.files("/standard/vol186/bergland-lab/Gio/subbaypass", "xtx", full.name=T)

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
  setkey(snps, MRK, subpool)
  setkey(xtx, MRK, subpool)
  xtx <- merge(xtx, snps)

  save(xtx, file="/scratch/aob2x/xtx_dest2.raw.Rdata")

### average
  xtx.ag <- xtx[,list(XtXst_median=median(XtXst), XtXst_mean=mean(XtXst),
                      neglogp_median=median(negLog10P), neglogp_mean=mean(negLog10P)), list(MRK, subpool)]
  xtx.ag[order(subpool)]

### SNPs
  snps <- fread("/standard/vol186/bergland-lab/Gio/subbaypass/dest_pos_table.txt")
  setnames(snps, "subpoolMRK", "MRK")

### merge
  setkey(snps, MRK, subpool)
  setkey(xtx.ag, MRK, subpool)
  xtx.ag <- merge(xtx.ag, snps)


##########
### contrast
##########

### get files
  fl <- list.files("/standard/vol186/bergland-lab/Gio/subbaypass", "contrast", full.name=T)

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
  setkey(snps, MRK, subpool)
  setkey(cont, MRK, subpool)
  xtx <- merge(cont, snps)

  save(cont, file="/scratch/aob2x/cont_dest2.raw.Rdata")

### average
  cont.ag <- cont[,list(C2_std_median=median(C2_std), C2_std_mean=mean(C2_std),
                        neglogp_median=median(negLog10P), neglogp_mean=mean(negLog10P)), list(MRK, subpool)]

  setkey(snps, MRK, subpool)
  setkey(cont.ag, MRK, subpool)
  cont.ag <- merge(cont.ag, snps)


  cont.ag[,p:=10^(-neglogp)]
  cont.ag[,q:=p.adjust(p, "fdr")]
  table(cont.ag$q<.005, cont.ag$chr, cont.ag$invName!="noInv")

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
  fisher.test(table(m$rnp<.05, m$XtXst>300))
  fisher.test(table(m$rnp<.05, m$q.y<.05))
  fisher.test(table(m$q.x<.0005, m$q.y<.05))

  fisher.test(table(m$q.x<.0005 & m$q.y<.05, m$q<.05))

  m[which.max(C2_std)]


###
