#### Collect and analyse the seasonality model

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(foreach)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions


log10_ceiling <- function(x) {
  10^(ceiling(log10(x)))
}


#####
#####
#####
#####
#####
### open metadata
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv")

samps <- fread("./dest_v2.samps_25Feb2023.csv")
### open filtering object
#keep.remove <- get(load("/project/berglandlab/DEST2.0_working_data/keep.fail.samps.Rdata"))
#pass.filt <- filter(keep.remove, keep == "PASS")$sampleId
#### open seasonal pairs
#seasonal.sets <- get(load("/project/berglandlab/DEST2.0_working_data/DEST2.seasonals.plusCore20.flip.met.Rdata"))
#set.samps <- filter(seasonal.sets, sampleId %in%  pass.filt)$sampleId

#### Central files
files.seas = system("ls ./GLM_omnibus_ALAN_APR122023", intern = T)
root = "/scratch/yey2sn/DEST2_analysis/seasonality/GLM_omnibus_ALAN_APR122023/"

#### binning analysis
seas.p.bin = 
  foreach(fil = files.seas, .combine = "rbind", .errorhandling = "remove" )%do%{
    
    message(fil)
    tmp <- get(load(paste(root, fil, sep = "/")))
    
    return(tmp)
    
  }

#### load object
snpdt.obj <- get(load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt.Rdata"))
setDT(snpdt.obj)

annotation <- get(load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_Cville_SNP_Metadata.Rdata"))
setDT(annotation)
names(annotation)[1:2] = c("chr","pos")
annotation$pos = as.numeric(annotation$pos)

left_join(seas.p.bin, snpdt.obj, by = c("chr", "pos")) %>%
  left_join(annotation, by = c("chr", "pos")) ->
  seas.p.bin.inv

#### #### #### 
#### 
#### 
#### 
#### Collect into permutation sets
#### Collect into permutation sets
#### Collect into permutation sets
#### Collect into permutation sets

outfolder = "/project/berglandlab/jcbnunez/Shared_w_Alan/GLM_omnibus_ALAN_APR122023_by_PERM"
### save as independent files
  foreach(p = 0:20)%do%{
    
    message(paste(p, 
                  sep = "/"))
    
    seas.p.bin %>% filter(perm == p 
    ) -> dat.p
    save(dat.p,
         file = paste(outfolder,
                      paste("GLM_out.perm_id", p, "Omni", "Rdata", sep = "."),
                      sep = "/"
                      ))
  }

#### Collect into permutation sets -- end
#### Collect into permutation sets -- end
#### Collect into permutation sets -- end
#### 
#### 
#### 
#### #### #### 
  

#### Analyze the p-value
seas.p.bin.inv -> dat.flt
ot_across_perms=
  foreach(mfet = "PhyloRan_LocRan", .combine = "rbind"  )%do%{
    
  foreach(p = 0:20, .combine = "rbind" )%do%{
    #foreach(chr = c("2L","2R","3L","3R"), .combine = "rbind")%do%{
    message(paste(p, 
                  #chr, 
                  sep = "/"))
    
    dat.flt %>% 
      filter(perm == p 
                       & model_features == mfet
    ) -> dat.low
    
    hist(dat.low$p_lrt, breaks = 100) -> hist.dat.low
    data.frame(hist.dat.low$mids ,    hist.dat.low$counts) %>%
      mutate(perm = p, 
             model_features = mfet
      )-> ot

    return(ot)
  } }

### Plot

setDT(ot_across_perms)

ot_across_perms %>%
  ggplot(
    aes(
      x=hist.dat.low.mids,
      y=(hist.dat.low.counts),
      color = perm == 0,
      group = perm
    )
  ) + 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +  geom_line(aes(alpha = perm == 0) ) +
  geom_vline(xintercept = 0.05) +
  facet_wrap(~model_features, scales = "free_y") + 
  theme_bw() +
  xlab("Seasonal P-value") +
  ylab("Count") +
  scale_alpha_manual(values = c(0.1, 1)) ->
  pvals.lines
#ggsave(pvals.lines, file = "pvals.lines.pdf", w = 9, h = 9)
ggsave(pvals.lines, file = "pvals.lines.png", w = 5, h = 4)

####
dat.flt %>% 
  group_by(perm, model_features, Putative_impact) %>%
  mutate(bin.id = case_when(
    p_lrt < 1e-5 ~ "e-5",
    p_lrt < 1e-4 ~ "e-4",
    p_lrt < 1e-3 ~ "e-3",
    p_lrt < 1e-2 ~ "e-2",
    p_lrt < 1e-1 ~ "e-1"
  )) ->
  dat.flt.bin.id

dat.flt.bin.id %>% 
  group_by(perm, bin.id, chr, model_features, Putative_impact) %>%
  summarize(N = n()) %>%
  .[complete.cases(.),] -> cum.dat
setDT(cum.dat)

  ggplot() +
    geom_boxplot(
      data = cum.dat[perm != 0][model_features == c("PhyloRan_LocRan")],
    aes(x=Putative_impact , color =chr ,y=log10(N)),
    outlier.shape = NA
    ) +
  geom_point(
    data = cum.dat[perm == 0][model_features == c("PhyloRan_LocRan")],
    aes(x=Putative_impact, color =chr ,y=log10(N)), size = 3) +
    coord_flip() +
    facet_wrap(chr~bin.id, scales = "free_x") ->
    cum.dis.plot
  ggsave(cum.dis.plot, file = "cum.dis.plot.pdf", w = 9, h  = 6)

####
#### ---> Window analysis
####
####
####
####
####
####
####  
  glm.out = seas.p.bin[model_features == "PhyloRan_LocRan"]
  
  glm.out %>%
    group_by(perm) %>%
    mutate(L = n()) %>% 
    mutate(rnp = rank(p_lrt)/L) ->
    glm.out
  
  glm.out %>%
    group_by(chr, pos, perm) %>%
    arrange(perm, pos) -> glm.out
  setDT(glm.out)
  
  glm.out.p = glm.out[perm == 0]
  
  win.bp <- 1e5
  step.bp <- 5e4
  
  setkey(glm.out.p, "chr")
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- glm.out %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)  
  
  ### start the summarization process
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     #win.i=1:10,
                     .errorhandling = "remove"#,
                     #.combine = "rbind"
  )%do%{
    
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.tmp <- glm.out %>%
                   filter(chr == wins[win.i]$chr) %>%
                   filter(pos >= wins[win.i]$start & pos <= wins[win.i]$end)
    setDT(win.tmp)
    
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]
    
    pr.i <- c(
      0.05
    )
    
    #tmpo <- foreach(
    # pr.i=thrs, 
    # .errorhandling="remove", 
    # .combine="rbind")%do%{
    win.tmp %>% 
      filter(!is.na(rnp), pr.i == pr.i ) %>%
      group_by(perm, chr) %>%
      summarise(pos_mean = mean(pos),
                pos_mean = mean(pos),
                pos_min = min(pos),
                pos_max = max(pos),
                win=win.i,
                pr=pr.i,
                rnp.pr=c(mean(rnp<=pr.i)),
                rnp.binom.p=c(binom.test(sum(rnp<=pr.i), 
                                         length(rnp), pr.i)$p.value),
                min.p.lrt=min(p_lrt),
                min.rnp=min(rnp),
                nSNPs = n(),
                sum.rnp=sum(rnp<=pr.i),
      )   %>%
      mutate(
        perm_type=ifelse(perm==0, "real","permuted"),
        invName=case_when(
          chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
          chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
          chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
          chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
          TRUE ~ "noInv"
        )) -> win.out
    
    win.out %>%
      group_by(perm_type,invName, chr, pos_mean ) %>%
      summarize(rnp.m05 = quantile(rnp.binom.p, 0.05)) ->
      out
    
    return(out)
    #}
    #tmpo
  }
  
  ##}
save(win.out, file = "win.out.Rdata")  
  
load("win.out.Rdata")  

win.out.dt = do.call("rbind", win.out)

win.out.dt %>%
  ggplot(
    aes(
      x=pos_mean/1e6,
      y=-log10(rnp.m05),
      color=perm_type
    )
  ) +
  geom_line() +
  facet_wrap(~chr, scales = "free_x", ncol = 1) ->
  rnp.plot 
ggsave(rnp.plot, file = "rnp.plot.pdf", w = 6, h = 7)

