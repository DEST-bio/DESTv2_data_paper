#### SNP-wise analyses -- PE Data

#### ==> Libs
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

#### ===> Load PE's data
PE_snps <- fread("erickson_top1.txt")

#### ===> Lift-over data
library(GenomicRanges)
library(rtracklayer)
chain_3to6="/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain"
chainObject <- import.chain(chain_3to6)

##iterate over positions
lift_over_top_snps=
foreach(i=1:dim(PE_snps)[1],
        .combine = "rbind")%do%{
          
          message(paste(i, dim(PE_snps)[1], sep = "|"))
          tmp = PE_snps[i]
          chr.i=paste("chr", tmp$chr, sep = "")
          pos.i=tmp$dm3.pos
          
          grObject <- GRanges(seqnames=c(chr.i), 
                              ranges=IRanges(start=pos.i, end=pos.i))
          results <- as.data.frame(liftOver(grObject, chainObject))
          
          data.frame(
            chr=tmp$chr,
            dm3.pos=tmp$dm3.pos,
            dm6.pos=results$start
          )
        }

save(lift_over_top_snps, file = "PE_outliers.lift_over_top_snps.Rdata")
# run liftOver

# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))

