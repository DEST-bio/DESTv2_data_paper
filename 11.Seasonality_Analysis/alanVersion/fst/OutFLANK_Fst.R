# ijob -A berglandlab -c20 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### capture command line arguments
  args = commandArgs(trailingOnly=TRUE)
  slurmJob=as.numeric(args[1])

### libraries
  #.libPaths(c("/scratch/aob2x/biol4559-R-packages-newer")); .libPaths()
  #install.packages("curl")
  library(curl)
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### source a few functions
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_input_freqs.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_correct.R")
  source("https://raw.githubusercontent.com/j-a-thia/genomalicious/master/R/outflank_mod_fst_nocorrect.R")
  # source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/OutFLANK.R")
  # source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Likelihood%20functions%20for%20OutFLANK.R")
  # source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")
  # source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/FST%20functions.R")

### load this function
  getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
    # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps[grepl("SRP002024", sampleId)]$sampleId

    ### filter to target
    seqSetFilter(genofile, variant.id=snps$id, sample.id=samples$sampleId)

    ### get annotations
    message("Annotations")
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snps$id, times=len1))

    # Extract data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                          list(variant.id=id)]

    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

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
    afi <- merge(af, snp.dt1.an, by="variant.id")
    afi <- merge(afi, snps, by.x="variant.id", by.y="id")

    afi[,af:=ad/dp]

    ### calculate effective read-depth
    afis <- merge(afi, samples[,c("sampleId", "nFlies", "locality",
                                  "lat", "long", "continent", "country", "province", "city",
                                  "min_day", "max_day", "min_month", "max_month", "year", "jday", "season",
                                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "subsample", "sampling_strategy",
                                  "SRA_Accession"), with=F], by="sampleId")

    afis[chr=="X|Y", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
    afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
    afis[,af_nEff:=round(af*nEff)/nEff]
    setnames(afis, "col", "annotation")
    ### return
    afis[,-c("n"), with=F]
  }

### open GDS file
  #genofile <- seqOpen("/scratch/aob2x/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.gds")
  genofile <- seqOpen("/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")
  genofile

### load meta-data file
  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")
  seasonal.sets <- samps[Recommendation=="Pass"][continent%in%c("Europe", "Asia")][set!="DrosRTEC"][,list(.N,
                                                                 season=c("spring", "fall"),
                                                                 sampleId=c(sampleId[which.min(jday)], sampleId[which.max(jday)]),
                                                                 diff=c(jday[which.min(jday)]-jday[which.max(jday)])),
                                                              list(country, locality, year, loc_rep, cluster2.0_k5)][N>1][abs(diff)>30]
  setkey(seasonal.sets, country, locality, year, loc_rep, cluster2.0_k5, sampleId)
  setkey(samps, country, locality, year, loc_rep, cluster2.0_k5, sampleId)
  seasonal.sets <- merge(seasonal.sets, samps)

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)
  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data] ### this is the global average frequency
  snps <- fread("/standard/vol186/bergland-lab/Gio/subbaypass/dest_pos_table.txt")
  setnames(snps, "subpoolMRK", "MRK")
  setkey(snps, chr, pos)
  snp.dt <- merge(snp.dt, snps[,c("variant.id", "invName"),with=F], by.y="variant.id", by.x="id")


### build a dictionary of windows
  win.bp <-  5e5
  step.bp <- 5e5

  setkey(snp.dt, "chr")

  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),.combine="rbind", .errorhandling="remove")%dopar%{
    # chr.i <- "2L"
    tmp <- snp.dt[J(chr.i)]
    data.table(chr=chr.i,
               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }

  wins[,window:=1:dim(wins)[1]]
  dim(wins)

### iterate thorugh windows to calculate a summary of your data
  setkey(snp.dt, chr, pos)

  # slurmJob <- 50
  out <- foreach(i=slurmJob, .errorhandling="remove")%do%{
    # i <- 120
    # message(i)
    ### get data for your sample for this window
      focalSNPs <-snp.dt[J(wins[i]$chr)][pos>=wins[i]$start & pos<=wins[i]$end]
      tmp.data <- getData(snps=focalSNPs, samples=seasonal.sets)
      tmp.data[,window:=i]
      tmp.data.clean <- tmp.data[!is.na(af_nEff)]

      tmp.data.clean[,locyear:=paste(locality, year, sep="_")]

      uniq_locyear <- unique(tmp.data.clean$locyear)
      fst.out <- foreach(u=uniq_locyear, .errorhandling="remove")%dopar%{
        # u <- uniq_locyear[1]
        snp.data <- tmp.data.clean[locyear==u]

        fst.out.tmp <- foreach(j=unique(snp.data$variant.id), .errorhandling="remove")%do%{
          # j<- unique(snp.data$variant.id)[36]
          message(paste(u, j, sep=" / "))
          #snp.data[variant.id==j]

          fst.cor <- outflank_mod_fst_correct(  snp.data[variant.id==j]$af_nEff,   snp.data[variant.id==j]$nEff, Ho=NULL)
          fst.nocor <- outflank_mod_fst_nocorrect(snp.data[variant.id==j]$af_nEff,   snp.data[variant.id==j]$nEff, Ho=NULL)

          mat <- matrix(c(snp.data[variant.id==j][order(season)]$af_nEff*snp.data[variant.id==j][order(season)]$nEff,
                          (1-snp.data[variant.id==j][order(season)]$af_nEff)*snp.data[variant.id==j][order(season)]$nEff),
                        nrow=2, ncol=2, byrow=T)

          fet <- fisher.test(mat)

          if(length(fst.cor)==1) {
            fstDF <- data.table(variant.id=unique( tmp.data$variant.id))
            return(fstDF, nPops=length(snp.data[variant.id==j]$af_nEff),
                  sp_allele_freq=snp.data[variant.id==j][season=="spring"]$af_nEff,
                  fa_allele_freq=snp.data[variant.id==j][season=="fall"]$af_nEff,
                  fet_p=fet$p.value, or=fet$estimate, lci=fet$conf.int[1], uci=fet$conf.int[2])
          } else {

            fstDF <- data.table(variant.id=j, nPops=length(snp.data[variant.id==j]$af_nEff), locyear=u,
                              cbind(as.data.frame(fst.nocor[c('FSTNoCorr', 'T1NoCorr', 'T2NoCorr')]), as.data.frame(fst.cor)),
                                    sp_allele_freq=snp.data[variant.id==j][season=="spring"]$af_nEff,
                                    fa_allele_freq=snp.data[variant.id==j][season=="fall"]$af_nEff,
                                    sp_allele_nEff=snp.data[variant.id==j][season=="spring"]$nEff,
                                    fa_allele_nEff=snp.data[variant.id==j][season=="fall"]$Eff,
                                    fet_p=fet$p.value, or=fet$estimate, lci=fet$conf.int[1], uci=fet$conf.int[2])
            return(fstDF)
          }
        }
        rbindlist(fst.out.tmp, fill=F)
      }
      fst.out <- rbindlist(fst.out, fill=F)


    ### format output
      out <- merge(fst.out,
            tmp.data[,list(annotation=unique(annotation), gene=unique(gene), meanRD=mean(af_nEff, na.rm=T), meanAF=mean(af_nEff, na.rm=T)), list(variant.id, chr, pos)],
            by="variant.id", all.y=T)
  }
  #out <- rbindlist(out, fill=T)

### save output
    outputdir="/scratch/aob2x/dest_fstoutput/"
    if(!dir.exists(outputdir)) {
      dir.create(outputdir)
    }
    save(out, file=paste(outputdir, "dest_fst_fet_window", slurmJob, ".Rdata", sep=""))
