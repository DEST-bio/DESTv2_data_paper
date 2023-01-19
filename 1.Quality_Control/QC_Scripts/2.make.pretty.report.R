### Pretty report Qulimap
### 

library(tidyverse)
library(magrittr)
library(vroom)
library(foreach)
library(patchwork)
library(reshape2)
## load dat
## 

#### ---> this is a file that contains the information of the files to fastQC, including the address.
qcsamps <-vroom("qcs.batch.txt", col_names = F)

#system(paste("ls", paste(qcsamps$X2[1], "raw_data_qualimapReport", sep = "/")),
#       intern = T) -> items.to.load

items.to.load <- c("coverage_across_reference.txt", "insert_size_across_reference.txt", "mapping_quality_across_reference.txt", "mapped_reads_gc-content_distribution.txt")

o.all <- foreach(k=1:length(items.to.load), .combine = "rbind")%do%{
tmp.o <- foreach(j=1:dim(qcsamps)[1], .combine = "rbind")%do%{
tmp <- vroom(paste(qcsamps$X2[j], "raw_data_qualimapReport", 
                               items.to.load[k], sep = "/")) %>%
  mutate(samp = qcsamps$X1[j],
         group = qcsamps$X3[j],
         metric = items.to.load[k])
names(tmp)[1:2] = c("pos","var")
return(tmp)
}

maxt = ifelse(items.to.load[k] %in% c("insert_size_histogram.txt",
                                      "duplication_rate_histogram.txt",
                                      "coverage_across_reference.txt"
                                      )
              , 1e9, 200)

tmp.o %>%
  filter(var != 0) %>%
  filter(var < maxt) %>%
  ggplot(aes(x=group, y=(var), fill=group)) + 
  ylab(paste(items.to.load[k])) + 
  geom_boxplot() +
  ggtitle("T-test: BAD vs REG", subtitle = t.test(var~ group,data = tmp.o)$p.value)-> tmp.box
  
tmp.o %>%
  filter(var > 0) %>%
  filter(var < maxt) %>%
  ggplot(aes(x=pos, y=(var), color=group)) + 
  ylab(paste(items.to.load[k])) + 
  geom_point() -> tmp.line

ggsave(tmp.box+tmp.line, file = paste(items.to.load[k],"png", sep ="."), h= 4, w=8)
return(tmp.o)
}

23513712+25286936+28110227+32079331+1348131+23542271 -> dmel.genome

o.all$metric = gsub(".txt", "", o.all$metric)

o.all %>% 
  filter(metric %in% c("coverage_across_reference", 
                       "insert_size_across_reference", 
                       "mapping_quality_across_reference")) %>%
  mutate(bam_slice = case_when(pos < dmel.genome ~ "Dmel",
                   pos >= dmel.genome ~ "Holo")) %>% 
  group_by(samp, bam_slice, metric) %>%
  summarise(mean = mean(var)) %>%
  reshape2::dcast(samp~metric+bam_slice, value.var = "mean") -> coverage.et.al
names(coverage.et.al)[1] = "Sample"

### Load the global mean

###  This is taken from MULTIQC!!
global <- vroom("All.bams.summary.txt")

### Merge
left_join(global, coverage.et.al) -> global.dmel.holo

#### 

qcs = qcsamps

dup.rat.head <- read.table( pipe( paste("sed -n 7p ", "/project/berglandlab/DEST/dest_mapped/DESTv2/" , 
                                        paste(qcs$X1[1], "/", sep = ""),
                                        paste(qcs$X1[1], ".mark_duplicates_report.txt", sep = ""), 
                                        sep = "")))

dup.rate <- foreach(i=1:dim(qcs)[1], .combine = "rbind")%do%{
  
  dup.rat <- read.table( pipe( paste("sed -n 8p ", "/project/berglandlab/DEST/dest_mapped/DESTv2/" , 
                                     paste(qcs$X1[i], "/", sep = ""),
                                     paste(qcs$X1[i], ".mark_duplicates_report.txt", sep = ""), 
                                     sep = ""))
  )
  
  dup.rat %<>% as.data.frame() %>% mutate(samp = qcs$X1[i])
  return(dup.rat)
}

names(dup.rate)[1:10] = dup.rat.head

names(qcs) = c("samp","path", "batch")

dup.rate %>%
  left_join(qcs) -> data.pcr.dups
####
###
data.pcr.dups %>%
  ggplot(aes(x=batch, y=PERCENT_DUPLICATION, fill = batch)) + 
  geom_boxplot() ->
  PERCENT_DUPLICATION
ggsave(PERCENT_DUPLICATION, file = "PERCENT_DUPLICATION.png", w = 5, h =4)

t.test(PERCENT_DUPLICATION~batch, data=data.pcr.dups)

#####

data.pcr.dups[,c("UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "samp")] -> data.pcr.dups.tmp
names(data.pcr.dups.tmp)[which(names(data.pcr.dups.tmp) == "samp")] = "Sample"

left_join(global.dmel.holo, data.pcr.dups.tmp) -> Mapping.stats.final

write.table(Mapping.stats.final, file = "Mapping.stats.final.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



