# ijob -A berglandlab -c5 -p largemem --mem=250G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.0.3; R

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(doMC)

### get allele frequencies
  getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
    # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps[grepl("SRP002024", sampleId)]$sampleId

    ### filter to target
    seqSetFilter(genofile, variant.id=snps$id, sample.id=samples$sampleId)

    ### get frequencies
    message("Allele Freqs")

    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")

    if(class(dp)[1]!="SeqVarDataList") {
      dp <- list(data=dp)
    }


    af <- data.table(ad=expand.grid(ad$data)[,1],
                     dp=expand.grid(dp$data)[,1],
                     sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                     variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
    message("merge")
    afi <- merge(af, snps, by.x="variant.id", by.y="id")

    afi[,af:=ad/dp]

    ### calculate effective read-depth
    afis <- merge(afi, samples[,c("sampleId", "set", "nFlies", "locality",
                                  "lat", "long", "continent", "country", "province", "city",
                                  "min_day", "max_day", "min_month", "max_month", "year", "jday",
                                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "subsample", "sampling_strategy",
                                  "SRA_Accession"), with=F], by="sampleId")

    afis[chr=="X|Y", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
    afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
    afis[,af_nEff:=round(af*nEff)/nEff]
    #setnames(afis, "col", "annotation")
    ### return
    afis[,-c("n"), with=F]
  }

### open GDS file
  #genofile <- seqOpen("/scratch/aob2x/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.gds")
  genofile <- seqOpen("/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")
  genofile

### load meta-data file
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency

### determine samples to use
  ### K8
    samps.use_k8 <- samps[,list(sampleId=sample(sampleId)), list(locality, k8=cluster2.0_k8, continent)]
    table(samps.use_k8$k8, samps.use_k8$continent)
    samps.use_k8 <- samps.use_k8[k8%in%c(2,7,8)][!continent%in%c("Africa", "Oceania")]
    setkey(samps.use_k8, locality, continent, sampleId)
    setkey(samps, locality, continent, sampleId)
    dim(samps.use_k8)
    samps.use_k8 <- merge(samps.use_k8, samps)
    dim(samps.use_k8)

### get data
  testVariant <- sort(sample(snp.dt[chr%in%c("2L", "2R", "3L", "3R")]$id, 5000))
  setkey(snp.dt, id)
  dat <- getData(snps=snp.dt[J(testVariant)], samples=samps.use_k8)
  dat <- merge(dat, samps.use_k8[,c("sampleId", "k8"), with=F], by="sampleId")
  dat[,list(nEff=mean(nEff, na.rm=T)), list(sampleId)][order(nEff)]

### weighted average
  wmean <- function(x, n) {
    sum(x*n, na.rm=T)/sum(n[!is.na(n)])
  }
  dat.ave <- dat[,list(af.neff.wmean=wmean(x=af_nEff, n=nEff),
                       neff.sum=sum(nEff, na.rm=T),
                        af.neff.mean=mean(af_nEff, na.rm=T),
                        neff.mean=mean(nEff, na.rm=T)),
                  list(chr, pos, variant.id, k8)]
  cor.test(dat.ave$af.neff.wmean, dat.ave$af.neff.mean)
  cor.test(dat.ave$af.neff.wmean, dat.ave$af.neff.mean)

  dat.ave.ag <- dat.ave[,list(nEff=mean(neff.mean)), list(k8)]

  
