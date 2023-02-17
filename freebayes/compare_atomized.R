### libraries
  library(data.table)
  library(SeqArray)

### convert Freebayes to GDS
  grep -v "#" ~/tst.atomized.vcf | cut -f1,2,3,4,5 > ~/sites.delim



### dl
  #system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/DEST/paramTest/snpcaller/dest.PoolSeq.PoolSNP.001.5.test.norep.ann.gds ~/.")
  #system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DESTv2_output/dest.PoolSeq.SNAPE.001.5.test.norep.ann.gds ~/.")


### open GDS files
  geno_ps <- seqOpen("~/dest.PoolSeq.SNAPE.001.5.test.norep.ann.gds")

### load freebayes snp table
  seqResetFilter(geno_ps)


  snp.ps <- data.table(chr=seqGetData(geno_ps, "chromosome"),
                        pos=seqGetData(geno_ps, "position"),
                        nAlleles=seqGetData(geno_ps, "$num_allele"),
                        ps_id=seqGetData(geno_ps, "variant.id"))
  snp.ps <- snp.ps[nAlleles==2]

  seqSetFilter(geno_ps, snp.ps$ps_id)
  snp.ps[,ps_af:=seqGetData(geno_ps, "annotation/info/AF")$data]

### atomized
  atom <- fread("~/sites.delim")
  setnames(atom, names(atom), c("chr", "pos", "class", "ref", "alt"))
  atom2 <- atom[!grep(",", alt)]

  atom[,nAlleles_fb:=str_count(alt, ",")+1]


### merge
  setkey(atom, chr, pos)
  setkey(snp.ps, chr, pos)
  m <- merge(atom, snp.ps, all=T)

  table(is.na(m$class), is.na(m$nAlleles))


### get data
  seqSetFilter(geno_fb, variant.id=snp.fb$fb_id[1])
  seqSetFilter(geno_ps, variant.id=4914)

  seqGetData(geno_fb, "annotation/format/DP")
  seqGetData(geno_fb, "annotation/format/AD")[[2]]

  seqGetData(geno_fb, "$ref")
  seqGetData(geno_fb, "$alt")
  seqGetData(geno_fb, "annotation/info/TYPE")

  table(seqGetData(geno_fb, "annotation/info/TYPE")$data)

  seqSetFilter(geno_fb, variant.id=snp.fb[pos==13613298]$fb_id)
