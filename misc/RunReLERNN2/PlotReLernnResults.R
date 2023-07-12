##Analysis of the estimated recombinates
args <- commandArgs(TRUE)
print(args[1])
relearn <- read.delim(args[1])
samplename <- sub(".*/(\\w+).*", "\\1", args[1])
print(samplename)
setwd("/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/results")
dir.create("/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/results/plots")
## Make subsets based on 'chrom' categories
subsets <- split(relearn, relearn$chrom)
for (subset in subsets) {
  plot(subset$start, subset$recombRate, type= "l", col="darkred", xlab = paste("Chromsosome", subset$chrom[1], "(Mb)"), main = paste(samplename), ylab="Recombination rate(c/bp)")
  png(filename = paste("/media/ssteindl/DESTv2_data_paper/misc/RunReLERNN2/results/plots/",samplename,"_", subset$chrom[1], ".png", sep = ""), width = 2000, height = 1200, units = "px", res = 300)
  }
dev.off()