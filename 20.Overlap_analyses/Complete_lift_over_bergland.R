#### Lift OVER bergland 
set = "bergland2016"
######
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)


#### ==> Libs
library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(ggvenn)
library(GenomicRanges)
library(rtracklayer)
chain_3to6="/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain"
chainObject <- import.chain(chain_3to6)

#####
load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2016.Bergland/6d_data.Rdata")

win.bp <- 1e5
step.bp <- win.bp+1

wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  tmp <- p %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
###963


####
wins[k] -> guide_ith

p %>%
  filter(chr == guide_ith$chr) %>%
  filter(pos >= guide_ith$start & pos < guide_ith$end) ->
  p_ith

#####
data = p_ith
#####

lift_over_pos=
  foreach(i=1:dim(data)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            
            message(paste(i, dim(data)[1], sep = "|"))
            tmp = data[i,]
            
            chr.i=paste("chr", tmp$chr, sep = "")
            pos.i=tmp$pos
            
            grObject <- GRanges(seqnames=c(chr.i), 
                                ranges=IRanges(start=pos.i, end=pos.i))
            results <- as.data.frame(liftOver(grObject, chainObject))
            
            data.frame(
              chr=tmp$chr,
              pos=tmp$pos,
              dm6.pos=results$start
            )
          }

lift_over_pos %>%
  mutate(SNP_id_dm6 = paste(chr,dm6.pos, sep = "_") ) %>%
  right_join(data) -> 
  lift_over_data

system(paste("mkdir ", set, "_LifOve", sep = "" ))

save(lift_over_data, 
     file = paste(paste(set, "_LifOve", sep = "" ),
                  "/",
                  paste(k,".",set,".Rdata", sep = ""), sep = ""
     ))

