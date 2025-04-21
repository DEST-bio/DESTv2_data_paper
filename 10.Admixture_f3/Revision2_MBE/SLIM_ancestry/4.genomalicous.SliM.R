#### Genomalicious on SLim Data
# module load Rgeospatial
library(genomalicious)
library(SeqArray)
library(foreach)
library(reshape2)
library(sp)

args = commandArgs(trailingOnly=TRUE)
job=as.numeric(args[1])

###bring in all simulation files
files2 <- system(paste("ls /gpfs2/scratch/jcnunez/DEST2.0_analysis/REVISION3_ADMIX/results/Agnostic.",job,".*", sep = ""),
                 intern = T)

o =
  foreach(i = c(files2),
          .combine = "rbind")%do%{
            
            replicate=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][2]
            pop=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][4]
            set=strsplit(strsplit(i,  "\\/")[[1]][8], "\\.")[[1]][1]
            
            tmp <- fread(i) %>% t()
            data.frame(AF=tmp) %>%
              mutate(pop=pop) %>%
              mutate(set=set) %>%
              mutate(replicate=replicate)
            
          }

########

for(k in 2:8){

o %>%
  filter(pop %in% c(0, 1, k)) ->
  tmp.f

names(tmp.f)[1:2] = c("FREQ","POOL")
# AF pop      set replicate
tmp.f %>% group_by(POOL) %>% 
  mutate(locus = 1:n()) %>%
  mutate(ref = "A", alt = "T") %>%
  mutate(ne = case_when(POOL == 0 ~ 25,
                        POOL == 1 ~ 25,
                        TRUE ~ 25)) %>%
  mutate(POOL = case_when(POOL == 0 ~ "A_parent1",
                   POOL == 1 ~ "B_parent2",
                   POOL == k ~ "C_derived",
                   ))->
  tmp.f

tmp.f %>% dcast(locus~POOL, value.var = "FREQ") %>%
  mutate(RMS = rowMeans(.[-1])) %>% 
  filter(RMS > 0.05) %>%
  .$locus -> flt.locus

tmp.f %>%
  filter(locus %in% flt.locus) ->
  tmp.f

setDT(tmp.f)
sfs_method="probs"

dadi <- dadi_inputs(
  tmp.f,
  type = "freqs",
  popCol = "POOL",
  locusCol = "locus",
  refCol = "ref",
  altCol = "alt",
  freqCol = "FREQ",
  indsCol = "ne",
  freqMethod = sfs_method
)

dadi <- na.omit(dadi) 

L=dim(dadi)[1]
####
### Write dadi file

fn <- paste("dadi_sim_objects/",
            sfs_method, ".",
            0, ".",
            1, ".",
            k, ".",
            "rep", ".", job, ".",
            "delim", sep="")

write.table(dadi,
            file=fn,
            sep="\t", quote=F, row.names=F)

### Metadata
fn2 <- paste("L_meta_sim_objects/",
             sfs_method, ".",
             0, ".",
             1, ".",
             k, ".",
             "rep", ".", job, ".",
             "meta", sep="")

an1_ne=25
an2_ne=25
Der_ne=25
write.table(data.frame(an1_ne, an2_ne, Der_ne, L),
            file=fn2,
            sep="\t", quote=F, row.names=F)
}
