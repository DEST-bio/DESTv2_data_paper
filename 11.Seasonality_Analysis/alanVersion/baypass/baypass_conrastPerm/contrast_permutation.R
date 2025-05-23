# ijob -A berglandlab_standard -c20 -p standard --mem=250G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### capture slurm id
  args <- commandArgs(trailingOnly=TRUE)
  subPool <- as.numeric(args[1]) #subPool=38

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
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
  setnames(snp.dt, "subpool", "subPool")

### c2_function
  c2_fun <- function(pij=pijs, cont, nperm=100) {
    # pij <- pijs; cont=cont.perm; nperm=100
    #create a nsnp x npop matrix with std allele frequencies
      message("standarize and unshrink")
      pij=t(matrix(pij$M_Pstd,nrow=max(pij$POP))) #be careful with transpose (important to remind that R always operates by column in matrix operations, see also below)

    #unshrink estimator
      pij=(pij-mean(pij))/sd(pij)

    ##create a matrix of (naive) permuted contrast
      #contrasts.mat=matrix(0,nrow=nperm+1,ncol=length(cont))
      #set.seed(101)

      #contrasts.mat[1,]=cont #first row = original contrast
      #message("randomize")
      #for(i in 2:nrow(contrasts.mat)){contrasts.mat[i,]=sample(cont)}
      ##contrasts.mat[2,]=cont #just pick up a random row inside contrasts.mat to check all goes ok in matrix operation below

    ### use pre-computed contrasts
      contrasts.mat <- as.matrix(dcast(data=cont, perm~sampleId, value.var="baypass.perm")[,-1])
      table(contrasts.mat[1,]==cont.perm[perm==0]$baypass.perm)

    ##compute all contrast all at once
      message("calculate c2")
      all.c2=((pij%*%t(contrasts.mat))**2)
      all.c2=t(t(all.c2)/rowSums(contrasts.mat**2)) #trick to speed computation as R operates by column in matrix operations

    ### format output
      message("expand grid")

      all.c2.dt <- as.data.table(all.c2)
      all.c2.dt[,MRK:=c(1:dim(all.c2)[1])]
      c2.dt <- melt(all.c2.dt, id.vars="MRK", variable.name="perm", value.name="C2_std")
      c2.dt[,perm:=as.numeric(gsub("V", "", perm))-1]


      c2.dt[,neglog10P:=-log10(pchisq(C2_std, 1, lower.tail=F, log=F))]

    ### return(
      return(c2.dt)

  }

##load contrast
  load(file="/standard/vol186/bergland-lab/alan/dest_baypass/contrast_pop.perm.Rdata")
  contrast=as.numeric(read.table("/standard/vol186/bergland-lab/Gio/dest2_season_contrast.txt")[1,])

### load standardized allele freq estimates

  fl <- list.files("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass", "yij", full.name=T)
  fl.use <- fl[grepl(paste("subpool_", subPool, "_", sep=""), fl)]

  c2.dt <- foreach(fl.i=fl.use)%do%{
    # fl.i <- fl[1]
    ### info
      file <- tstrsplit(fl.i, "/")%>%last()
      subPool <- tstrsplit(file, "_")[[2]]%>%as.numeric()
      rep <- tstrsplit(file, "_")[[3]]%>%as.numeric()

    ### load
      message(file)
      pijs=fread(fl.i, data.table = F)

    ### run function
      c2.dt <- c2_fun(pij=pijs, cont=cont.perm)

      c2.dt[,subPool:=subPool]
      c2.dt[,rep:=rep]

    ### return
      return(c2.dt)
  }
  c2.dt <- rbindlist(c2.dt)
  c2.dt.ag <- c2.dt[,list(C2_std_median=median(C2_std), neglog10P_median=median(neglog10P), .N), list(perm, MRK, subPool)]
  c2.dt.ag
  setkey(c2.dt.ag, subPool, MRK)
  setkey(snp.dt, subPool, MRK)
  c2.dt.ag <- merge(c2.dt.ag, snp.dt)

  summary(c2.dt.ag[,list(pr=mean(C2_std_median[perm==0] > C2_std_median[perm!=0] )), list(chr, pos)])

  c2.dt.ag[,list(thr=quantile(C2_std_median[perm!=0], .995))]
  prop.table(table(c2.dt.ag[perm==0]$C2_std_median> 6.04))

### save()
  save(c2.dt.ag, file=paste("/scratch/aob2x/dest2_baypass/contrast_perms_v3/contrast_perm_subpool_", subPool, ".Rdata", sep=""))
  save(c2.dt,    file=paste("/scratch/aob2x/dest2_baypass/contrast_perms_v3/RAW_contrast_perm_subpool_", subPool, ".Rdata", sep=""))
