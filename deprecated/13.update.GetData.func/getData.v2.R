#### Get Data program ####

getData <- function(chr="2L", start=14617051, end=14617051, samplesUse=samps$sampleId) {
  # chr="2L"; start=14617051; end=14617051; samplesUse=samps$sampleId
  
  ### filter to target
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id, sample.id=samplesUse, verbose=T)
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  if(class(dp)[1]!="SeqVarDataList") {
    dp.list <- list()
    dp.list$data <- dp
    dp <- dp.list
  }
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp$data)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  #afi <- merge(af, snp.dt1.an, by="variant.id")
  afi <- merge(af, snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  #afis[,-c("n"), with=F]
  afis[,c("sampleId", "af_nEff", "nEff"), with=F]
  
}
