### libraries
  library(data.table)
  library(SeqArray)

### convert Freebayes to GDS
  vcf.fn="/Users/alanbergland/bam.tst.2L_13600000-13700000.freebayes.vcf"
  gds.fn=gsub(".vcf", ".gds", vcf.fn)

  seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", verbose=T, optimize=T)

### dl
  system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.norep.ann.gds ~/.")


### load freebayes snp table
  snp.dt <- 
