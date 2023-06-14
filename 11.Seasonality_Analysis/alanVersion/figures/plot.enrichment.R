# ijob -A berglandlab_standard -c5 -p standard --mem=4G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)

### compile data
  fl <- list.files("/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/enrichment", full.names=T)

  oo <- foreach(fl.i=fl, .combine="rbind")%do%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    oo
  }

  save(oo, file="~/mega_enrichment.Rdata")

### plot
  library(ggplot2)
  library(data.table)

  system("scp aob2x@rivanna.hpc.virginia.edu:~/mega_enrichment.Rdata ~/.")
  load("~/mega_enrichment.Rdata")

### collapse

  oo <- oo[,list(N=sum(N)), list(min_p, max_p, perm, pops, model_features, chr)]

  oo[,model_features:=factor(model_features, levels=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))]

  oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, max_p, model_features, pops)]

  oo <- rbind(oo, oo.ag, fill=T)

####
  foreach(pops.i=unique(oo$pops))%do%{
    message(pops.i)

    p1 <-
    ggplot(data=oo[perm!=0][chr=="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
    geom_line(data=oo[perm==0][chr=="genome"][pops==pops.i], aes(x=log10(max_p), y=(N), group=perm), color="red") +
    facet_grid(model_features~chr, scales="free_y")

    p2 <-
    ggplot(data=oo[perm!=0][chr!="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
    geom_line(data=oo[perm==0][chr!="genome"][pops==pops.i], aes(x=log10(max_p), y=(N), group=perm), color="red") +
    facet_grid(model_features~chr, scales="free_y")

    layout <-
    "ABBBB"

    mega <- p1 + p2 + plot_layout(design=layout) +  plot_annotation(title = pops.i)
    ggsave(mega, file=paste("~/modelEnrichment.new.", pops.i, ".png", sep=""), h=8, w=8)

  }
