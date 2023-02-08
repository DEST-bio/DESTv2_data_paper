### libraries
  library(data.table)
  library(SeqArray)

### convert Freebayes to GDS
  #vcf.fn="/Users/alanbergland/bam.tst.2L_13600000-13700000.freebayes.vcf"
  #gds.fn=gsub(".vcf", ".gds", vcf.fn)
#
  #seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", verbose=T, optimize=T)

### dl
  #system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/DEST/paramTest/snpcaller/dest.PoolSeq.PoolSNP.001.5.test.norep.ann.gds ~/.")

### open GDS files
  geno_fb <- seqOpen("/Users/alanbergland/bam.tst.2L_13600000-13700000.freebayes.gds")
  geno_ps <- seqOpen("~/dest.PoolSeq.PoolSNP.001.5.test.norep.ann.gds")

### load freebayes snp table
  seqResetFilter(geno_fb)
  seqResetFilter(geno_ps)
  snp.fb <- data.table(chr=seqGetData(geno_fb, "chromosome"),
                        pos=seqGetData(geno_fb, "position"),
                        nAlleles=seqGetData(geno_fb, "$num_allele"),
                        fb_id=seqGetData(geno_fb, "variant.id"))

  snp.ps <- data.table(chr=seqGetData(geno_ps, "chromosome"),
                        pos=seqGetData(geno_ps, "position"),
                        nAlleles=seqGetData(geno_ps, "$num_allele"),
                        ps_id=seqGetData(geno_ps, "variant.id"))
  snp.ps <- snp.ps[nAlleles==2]

  seqSetFilter(geno_ps, snp.ps$ps_id)
  snp.ps[,ps_af:=seqGetData(geno_ps, "annotation/info/AF")$data]

### merge
  setkey(snp.fb, chr, pos)
  setkey(snp.ps, chr, pos)
  m <- merge(snp.fb, snp.ps, all=T)

  table(m$nAlleles.x, m$nAlleles.y)


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
