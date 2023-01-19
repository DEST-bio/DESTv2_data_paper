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

  save(o, file="~/DESTv2_data_paper/first_batch_qc/depth_per_contig.Rdata")
