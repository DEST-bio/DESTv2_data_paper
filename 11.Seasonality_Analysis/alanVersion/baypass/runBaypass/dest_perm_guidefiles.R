# ijob -A berglandlab -c1 -p standard --mem=10G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(data.table)
  library(foreach)

### sample metadata
  samps <- fread("/standard/vol186/bergland-lab/DEST_v2/seasonal_sets_DESTv2.csv")
  samps2 <- fread("/standard/vol186/bergland-lab/Gio/dest_v2.samps_8Jun2023.csv")
  samps <- merge(samps, samps2, by="sampleId")

### make the permuted contrast file
  load(file="/standard/vol186/bergland-lab/alan/dest_baypass/contrast_pop.Rdata")
  cont <- merge(cont, samps[,c("sampleId", "locality.x", "year.x"), with=F], by="sampleId")

### permutaion within locality:year
  set.seed(1234)
  cont.perm <- foreach(i=1:1000, .combine="rbind")%do%{
    cont[,list(real=old, perm=sample(old, replace=F), i=i), list(locality.x, year.x)]
  }

  cont.perm.ag <- cont.perm[,list(diff=sum(real!=perm), .N), list(i)]
  cont.perm.ag[,frac:=diff/N]
  cont.perm.ag[order(diff)]
  cont.perm.ag[,rank:=rank(-frac, ties="first")]
  setkey(cont.perm, i)

  foreach(x=1:10)%do%{
    # x<- 1
    tmp <- cont.perm[J(cont.perm.ag[rank==x]$i)]
    table(tmp$real==tmp$perm)
    write.table(tmp$perm, quote=F, row.names=F, sep="\t", col.names=F, file=paste("/standard/vol186/bergland-lab/Gio/dest2_season_contrast_", x, ".txt", sep=""))
    paste("/standard/vol186/bergland-lab/Gio/dest2_season_contrast_", x, ".txt", sep="")
  }

### modify the control file to include permutation replicates
  guide <- fread("/standard/vol186/bergland-lab/alan/dest_baypass/subpoolcontrol.txt", header=F)

  guide.perm <- guide[,list(perm=1:10), list(V1, V2)]

  write.table(guide.perm, quote=F, row.names=F, col.names=F, sep=" ", file="/standard/vol186/bergland-lab/alan/dest_baypass/subpoolcontrol_perm.txt")
  guide
