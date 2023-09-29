### libraries
  library(data.table)
  library(ggplot2)
  library(SeqArray)

### load sample meta-data
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_8Jun2023.csv")

### open GDS file
  genofile <- seqOpen("/Users/alanbergland/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")

### seasonal set
  load("~/seasonalpair.pca.meta.Rdata")
  ss_noCore20 <- seasonalpair.pca.meta[delta.T.mag=="Steep"][delta.T.sign==-1][Core20_sat==F]$sampleId
  ss_Core20 <- seasonalpair.pca.meta[Core20_sat==T]$sampleId

### function
  getData <- function(variant, samps2use=seasonal.sets$sampleId) {
    # variant=snp.dt[chr=="2L"][pos==14617051]$variant.id;
    # samps2use=seasonal.sets$sampleId
    # variant=595225

    ### filter to target
      seqResetFilter(genofile)
      seqSetFilter(genofile, variant.id=variant, sample.id=samps2use)

    ### get frequencies
      message("Allele Freqs")

      ad <- seqGetData(genofile, "annotation/format/AD")$data
      dp <- seqGetData(genofile, "annotation/format/DP")

      af <- data.table(ad=expand.grid(ad)[,1],
                       dp=expand.grid(dp)[,1],
                       sampleId=rep(seqGetData(genofile, "sample.id"),
                                    dim(ad)[2]),
                       variant.id=rep(seqGetData(genofile, "variant.id"),
                                      each=dim(ad)[1]),
                        chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"))

    ### tack them together
      af[,af:=ad/dp]

    ### calculate effective read-depth
      afis <- merge(af, samps[,c("sampleId", "nFlies"), with=F], by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
      afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      afis
  }

### load in model output
  load("/Users/alanbergland/NoCore20_NewCore20_OrigCore20.Rdata")

### evaluate
  dat <- getData(variant=mm2[perm==0][which.min(p_lrt.y)]$variant.id.x[1], samps2use=ss_noCore20)
  dat <- merge(dat, seasonalpair.pca.meta, by="sampleId")

  ggplot(data=dat, aes(x=season, y=af_nEff, group=loc.y, color=locality)) + geom_line()
