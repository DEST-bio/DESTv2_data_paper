#ijob -c20 --mem=20G -p standard -A berglandlab

# module load gcc/9.2.0 openmpi/3.1.6 R/4.2.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### load output
  fl <- list.files("/scratch/aob2x/depth/", "summary.txt")
  o <- foreach(fl.i=fl, .combine="rbind")%do%{
    # fl.i <- fl[1]
    dt <- fread(paste("/scratch/aob2x/depth/", fl.i, sep=""))
    dt[,samp:=gsub(".mosdepth.summary.txt", "", fl.i)]
  }

  setkey(o, chrom)
  o <- o[J(c(c("2L", "2R", "3L", "3R"), paste("sim", c("2L", "2R", "3L", "3R"), sep="_")))]

  save(o, file="~/depth_per_contig.Rdata")
